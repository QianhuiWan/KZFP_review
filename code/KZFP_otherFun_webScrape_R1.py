import requests
import csv
import re
import os
import argparse
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm

# we can customize keywords by ourselves
repeat_keywords = [
    "repeat", "transposon", "retrotransposon", "ERV", "LINE", "SINE",
    "TE", "mobile element", "heterochromatin"
]

other_keywords = [
    "immune", "immunity", "inflammation", "development", "differentiation",
    "apoptosis", "cell cycle", "pluripotency", "neural", "signaling",
    "pathway", "tumor", "cancer", "stress", "stem cell", "metabolism"
]

def clean_text(text):
    return re.sub(r'\s+', ' ', text.strip().lower())

def summarize(abstract):
    sentences = abstract.split('. ')
    return '. '.join(sentences[:2]) + '.' if len(sentences) > 1 else abstract

def classify_focus_with_score(abstract):
    repeat_score = sum(abstract.count(kw.lower()) for kw in repeat_keywords)
    other_score = sum(abstract.count(kw.lower()) for kw in other_keywords)

    if repeat_score > 1.5 * other_score:
        focus = "Repeat regulation"
    elif other_score > 1.5 * repeat_score:
        focus = "Other function"
    else:
        focus = "Both / Unclear"

    review_needed = "Yes" if abs(repeat_score - other_score) < 2 else "No"
    return focus, repeat_score, other_score, review_needed

def build_url(paper):
    if 'pmcid' in paper:
        return f"https://europepmc.org/article/PMC/{paper['pmcid']}"
    elif 'id' in paper and paper.get('source') == 'MED':
        return f"https://europepmc.org/article/MED/{paper['id']}"
    else:
        return "N/A"

def fetch_all_kzfp_papers():
    print("Fetching all KZFP-related function papers from Europe PMC...")
    query = '("KZFP" OR "KRAB zinc finger protein" OR "KRAB-ZFP" OR "KRAB-ZNF" OR "KRAB ZNF" OR "KRAB-ZNF protein") AND (function* OR role* OR pathway* OR interact* OR express* OR regulat* OR mechanis* OR immune* OR cancer* OR disease* OR develop*)'
    base_url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
    page_size = 100
    cursor = "*"
    papers = []
    seen_ids = set()

    with tqdm(desc="Downloading pages") as pbar:
        while True:
            params = {
                'query': query,
                'format': 'json',
                'pageSize': page_size,
                'cursorMark': cursor,
                'resultType': 'core',
                'sort': 'cited desc'
            }

            response = requests.get(base_url, params=params)
            data = response.json()
            results = data.get("resultList", {}).get("result", [])
            cursor = data.get("nextCursorMark", None)

            if not results or cursor is None:
                break

            for r in results:
                uid = r.get('id')
                if uid in seen_ids:
                    continue
                seen_ids.add(uid)

                title = r.get('title', 'N/A')
                authors = r.get('authorString', 'N/A')
                journal = r.get('journalTitle', 'N/A')
                year = r.get('pubYear', 'N/A')
                doi = r.get('doi', 'N/A')
                abstract_raw = r.get('abstractText', '')
                abstract = clean_text(abstract_raw)
                summary = summarize(abstract_raw)
                focus, repeat_score, other_score, review_needed = classify_focus_with_score(abstract)
                url = build_url(r)

                papers.append([
                    title, authors, year, journal, doi, summary,
                    focus, repeat_score, other_score, review_needed, url
                ])

            pbar.update(1)

    return papers

def save_to_tsv(papers, output_dir):
    output_path = os.path.join(output_dir, 'kzfp_all_papers.tsv')
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'Title', 'Authors', 'Year', 'Journal', 'DOI', 'Summary',
            'Predicted Focus', 'Repeat Score', 'Other Score', 'Review Needed', 'EuropePMC URL'
        ])
        for paper in papers:
            writer.writerow(paper)

def plot_focus_distribution(papers, output_dir):
    focus_counts = Counter([p[6] for p in papers])  # column 6 = Predicted Focus
    labels = list(focus_counts.keys())
    values = list(focus_counts.values())

    plt.figure(figsize=(8, 5))
    plt.bar(labels, values)
    plt.title("Predicted Functional Focus of KZFP Papers")
    plt.ylabel("Number of Papers")
    plt.xticks(rotation=20)
    plt.tight_layout()

    plot_path = os.path.join(output_dir, "kzfp_focus_distribution.pdf")
    plt.savefig(plot_path)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Classify KZFP papers from Europe PMC and save results.")
    parser.add_argument('--output_dir', type=str, default='kzfp_output',
                        help='Directory to save output files (default: kzfp_output)')
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    papers = fetch_all_kzfp_papers()
    print(f"Total papers fetched and classified: {len(papers)}")

    save_to_tsv(papers, output_dir)
    print(f"Saved full results to '{output_dir}/kzfp_all_papers.tsv'")

    plot_focus_distribution(papers, output_dir)
    print(f"Bar plot saved to '{output_dir}/kzfp_focus_distribution.pdf'")



