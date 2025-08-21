import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

# ==== config ====
CYTOBAND_FILE = "cytoBand.txt"
KZFP_BED = "KZFPs.bed"
CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]
LABEL_STYLE = "chrband"   # 可选: chrband, genomic, both
OUTPUT = "KZFP_karyotype.svg"

# ==== 读取 cytoband 数据 ====
cyto = pd.read_csv(
    CYTOBAND_FILE,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "band", "stain"]
)
cyto = cyto[cyto["chrom"].str.replace("chr", "").isin(CHROM_ORDER)]
cyto["chrom"] = cyto["chrom"].str.replace("chr", "")

# 每条染色体长度
chrom_sizes = cyto.groupby("chrom")["end"].max()

# 分染色体聚合 cytoband 数据
cyto_sorted = cyto.groupby("chrom").apply(lambda df: df.sort_values("start")).reset_index(drop=True)

# ==== 读取 KZFP BED 文件 ====
kzfp = pd.read_csv(KZFP_BED, sep="\t", header=None,
                   names=["chrom", "start", "end", "name"])
kzfp["chrom"] = kzfp["chrom"].str.replace("chr", "")

# ==== 绘图 ====
fig, ax = plt.subplots(figsize=(16, 12))
ax.set_xlim(0, 15)
ax.set_ylim(0, 25)
ax.axis("off")

# 定义每列最多几个染色体
cols = 6
rows = (len(CHROM_ORDER) + cols - 1) // cols

# 各种染色体条带颜色
color_map = {
    "gneg": "white", "gpos25": "#C0C0C0", "gpos50": "#888888",
    "gpos75": "#444444", "gpos100": "black",
    "acen": "#FF0000", "gvar": "#E0E0E0", "stalk": "#CCCCCC"
}

# 绘制每条染色体
chrom_positions = {}
for idx, chrom in enumerate(CHROM_ORDER):
    col = idx // rows
    row = idx % rows
    x0 = col * 2.2
    y0 = 23 - row * 1.0

    chrom_positions[chrom] = (x0, y0)

    chrom_cyto = cyto_sorted[cyto_sorted["chrom"] == chrom]
    chrom_len = chrom_sizes[chrom]

    for _, r in chrom_cyto.iterrows():
        start = r["start"] / chrom_len
        end = r["end"] / chrom_len
        rect = patches.Rectangle(
            (x0, y0 + start), 0.6, end - start,
            facecolor=color_map.get(r["stain"], "white"),
            edgecolor="black", lw=0.2
        )
        ax.add_patch(rect)

    ax.text(x0 + 0.8, y0 + 0.5, chrom, fontsize=10, weight="bold")

# ==== 绘制 KZFP 标签 ====
for _, row in kzfp.iterrows():
    chrom = row["chrom"]
    if chrom not in chrom_positions:
        continue
    pos = (row["start"] + row["end"]) / 2
    chrom_len = chrom_sizes[chrom]
    frac = pos / chrom_len
    x0, y0 = chrom_positions[chrom]
    label_y = y0 + frac
    label_x = x0 + 0.7

    # 标签文字
    label = row["name"] if LABEL_STYLE == "gene" else f"{chrom}{row['name']}"
    ax.plot([x0+0.6, label_x], [label_y, label_y], color="red", lw=0.5)
    ax.text(label_x + 0.1, label_y, label, fontsize=6, color="red", va="center")

plt.tight_layout()
plt.savefig(OUTPUT, dpi=300)
plt.show()

print(f"✅ KZFP 核型图已生成: {OUTPUT}")


