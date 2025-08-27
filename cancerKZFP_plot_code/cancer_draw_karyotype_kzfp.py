# "get karyotype_KZFP plot"

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from bisect import bisect_left

# ===== config =====
## UCSC cytoBand (need to download first)
CYTO_BED = "~/githubRepo/KZFP_review/input_data/hg38_cytoBand.txt"   
## same genome version as cytoband（e.g.hg38）   
CHROM_SIZES = "~/githubRepo/KZFP_review/input_data/filtered_hg38_p14.chrom.sizes"   
## KZFP location at hg38 genome
KZFP_BED = "~/githubRepo/KZFP_review/input_data/kzfps_hg38_cancer_geneCoordinates.bed"     
## output path + output file name    
OUTPUT = "~/githubRepo/KZFP_review/outputs/cancer_KZFPs.pdf"
## 'cytoband' | 'genomic' | 'both'
LABEL_MODE = "cytoband"        

# chr order
CHROMS = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]

# ===== read in data =====
cyto = pd.read_csv(
    CYTO_BED, sep=r"\s+", header=None,
    names=["chrom","start","end","band","stain"]
)
cyto = cyto[cyto["chrom"].isin(CHROMS)].copy()

sizes = pd.read_csv(
    CHROM_SIZES, sep=r"\s+", header=None, names=["chrom","size"]
)
sizes = sizes[sizes["chrom"].isin(CHROMS)]
chrom_size = dict(sizes.values)


# create cytoband index
cyto_idx = {}
for c, sub in cyto.groupby("chrom", sort=False):
    sub = sub.sort_values("start").reset_index(drop=True)
    cyto_idx[c] = {
        "starts": sub["start"].tolist(),
        "ends":   sub["end"].tolist(),
        "bands":  sub["band"].tolist()
    }


def pos_to_cytoband(chrom, pos):
    if chrom not in cyto_idx: return None
    starts = cyto_idx[chrom]["starts"]
    i = bisect_left(starts, pos)
    idx = i-1 if i>0 else 0
    if cyto_idx[chrom]["starts"][idx] <= pos < cyto_idx[chrom]["ends"][idx]:
        return cyto_idx[chrom]["bands"][idx]
    return None

def label_for(chrom, pos):
    if LABEL_MODE == "genomic":
        return f"{chrom}:{pos/1e6:.1f}Mb"
    band = pos_to_cytoband(chrom, pos)
    cyt = f"{chrom}{band}" if band else chrom
    if LABEL_MODE == "cytoband":
        return cyt
    return f"{cyt} ({chrom}:{pos/1e6:.1f}Mb)"

# read in KZFP
kzfp = pd.read_csv(
    KZFP_BED, sep=r"\s+", header=None,
    names=["chrom","start","end","name"]
)
kzfp = kzfp[kzfp["chrom"].isin(CHROMS)].copy()
kzfp["mid"] = ((kzfp["start"] + kzfp["end"])//2).astype(int)

# ===== plotting（simple version）=====
# parameters（column * row）
cols = 6
rows = (len(CHROMS) + cols - 1)//cols
fig_w, fig_h = 16, 12
plt.figure(figsize=(fig_w, fig_h))
ax = plt.gca()
ax.axis("off")

# chromosome display
W = 0.55          # chr width
XGAP = 2.2        # column gaps
YGAP = 1.0        # row gaps
top = rows*YGAP + 1.0

# save chr axis position on plot for later use 
chrom_box = {}

for i, chrom in enumerate(CHROMS):
    col = i // rows
    row = i % rows
    x0 = col*XGAP
    y0 = top - row*YGAP - 0.1  # white space up
    chrom_box[chrom] = (x0, y0)

    # chr length 
    L = chrom_size[chrom]

    # band to vertical lines
    sub = cyto[cyto["chrom"]==chrom].sort_values("start")
    for _, r in sub.iterrows():
        s, e = r["start"], r["end"]
        frac_s, frac_e = s/L, e/L
        # unit=1, chr length
        rect = patches.Rectangle((x0, y0+frac_s), W, (frac_e-frac_s),
                                 fill=False, linewidth=0.6)
        ax.add_patch(rect)

    # chr number
    ax.text(x0 + W + 0.15, y0 + 0.5, chrom.replace("chr",""),
            fontsize=9, va="center")

# ===== label KZFP =====
# avoid overlapping
from collections import defaultdict
y_taken = defaultdict(list)  # chrom -> used y positions

for _, g in kzfp.sort_values(["chrom","mid","name"]).iterrows():
    chrom, pos, name = g["chrom"], int(g["mid"]), g["name"]
    if chrom not in chrom_box or chrom not in chrom_size: 
        continue
    L = chrom_size[chrom]
    frac = max(0, min(1, pos / L))
    x0, y0 = chrom_box[chrom]
    y = y0 + frac
    x_tick = x0 + W

    # 文本与引线位置（简单避让：若同一 y 已占用，向上/下微移）
    dy = 0.0
    step = 0.02
    while any(abs((y+dy) - yy) < 0.02 for yy in y_taken[chrom]):
        dy += step if (len(y_taken[chrom])%2==0) else -step
    y_final = y + dy
    y_taken[chrom].append(y_final)

    # 引线（不设颜色，采用默认）
    ax.plot([x_tick, x_tick+0.6], [y, y_final], linewidth=0.6)
    # 标签
    txt = f"{name} {label_for(chrom, pos)}"
    ax.text(x_tick+0.65, y_final, txt, fontsize=6, va="center")

plt.tight_layout()
plt.savefig(OUTPUT, dpi=300)
print(f"Saved: {OUTPUT}")




