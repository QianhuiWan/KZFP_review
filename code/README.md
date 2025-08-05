# KZFP Literature Classifier

This Python script fetches and classifies research papers related to **KRAB zinc finger proteins (KZFPs)** from Europe PMC. It categorizes each paper based on its main research focus (e.g., repeat regulation vs. other biological functions), and outputs a clean, tabular result for downstream analysis.

---

## What it does

- Queries Europe PMC for KZFP-related function papers (e.g., regulation, immunity, development).
- Parses abstracts and classifies them as:
  - **Repeat regulation**
  - **Other function**
  - **Both / Unclear**
- Outputs:
  - A TSV file with metadata, summary, classification, keyword match scores.
  - A PDF bar plot showing counts for each focus category.
  - A flag indicating whether manual review is needed (if classification confidence is low).
  - Direct links to abstracts on [EuropePMC](https://europepmc.org/).

---

## Requirements

Install required Python packages:

```bash
pip install requests matplotlib tqdm
```

---

## Usage

```bash
python3 ~/Users/qwan/githubRepo/coh_bioLLM/tests_QW_R2/KZFP_func/KZFP_otherFun_webScrape.py --output_dir your/path/to/output_results_kzfp
```

- `--output_dir`: Folder to save output files. Defaults to `./kzfp_output`.

---

## Output files

| File | Description |
|------|-------------|
| `kzfp_all_papers.tsv` | All classified papers in tab-separated format |
| `kzfp_focus_distribution.pdf` | Bar chart summarizing paper counts by predicted focus |

Each row in `kzfp_all_papers.tsv` includes:

- Title
- Authors
- Year, Journal, DOI
- Abstract Summary
- Predicted Focus (`Repeat regulation`, `Other function`, `Both / Unclear`)
- Keyword Match Scores (`Repeat Score`, `Other Score`)
- `Review Needed` flag
- EuropePMC URL

---

## Customization

### Keyword Matching

You can modify the keyword lists in the script:

```python
repeat_keywords = [...]
other_keywords = [...]
```

These are used to assign focus categories based on abstract content.

---

## Example Plot

*The PDF plot `kzfp_focus_distribution.pdf` will be saved in your output folder.*

---

## Example Use Cases

- Analyze trends in KZFP research topics
- Identify underexplored biological functions (what we are interested)
- Prioritize papers for further manual review or full-text reading

---

## Limitations

- Classification is based only on abstract text — some papers may require manual review for accuracy.
- EuropePMC API may rate-limit very large queries; this script paginates responsibly.

---

## Acknowledgments

Built using the [Europe PMC RESTful API](https://europepmc.org/RestfulWebService).

---

## Contact

For feedback or issues, please open an issue or contact Qianhui Wan.


