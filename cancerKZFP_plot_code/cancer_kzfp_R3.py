#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, sys, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(description="Plot genome ideogram from cytoBand.txt and optional BED labels")
    ap.add_argument("--cyto", required=True, help="UCSC cytoBand.txt (gz ok)")
    ap.add_argument("--bed", default=None, help="Optional BED: chrom start end name")
    ap.add_argument("--genome", default="hg38", help="Just for title, e.g., hg38/mm10")
    ap.add_argument("--out", default="ideogram.png", help="Output image")
    ap.add_argument("--cols", type=int, default=6, help="Columns of layout")
    return ap.parse_args()

def load_cyto(path):
    # UCSC cytoBand.txt.gz columns:
    # chrom, chromStart, chromEnd, name, gieStain
    cyto = pd.read_csv(path, sep="\t", header=None,
                       names=["chrom","chromStart","chromEnd","name","gieStain"])
    cyto = cyto.rename(columns={"chromStart":"start","chromEnd":"end"})
    return cyto

def human_chrom_order(chroms):
    # keep typical human order
    wanted = [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]
    # filter by existence
    ordered = [c for c in wanted if c in set(chroms)]
    # plus any remaining chromosomes at the end (e.g., chrM)
    rest = [c for c in chroms if c not in set(ordered)]
    return ordered + rest

def band_color(g):
    # UCSC gieStain -> facecolor
    # gneg: white; gpos25..100: grayscale; acen: centromere (red); gvar/stalk: light gray
    if g.startswith("gpos"):
        level = int(g.replace("gpos",""))
        gray = 1.0 - level/100.0
        return (gray, gray, gray)
    if g == "gneg":
        return (1,1,1)
    if g in ("gvar","stalk"):
        return (0.85,0.85,0.85)
    if g == "acen":
        return (0.80,0.25,0.25)  # reddish
    # fallback
    return (0.9,0.9,0.9)

def label_for_band(cyto_df, chrom, pos):
    sub = cyto_df[cyto_df["chrom"]==chrom]
    hit = sub[(sub["start"]<=pos) & (pos<sub["end"])]
    if not hit.empty:
        return hit.iloc[0]["name"]
    return ""

def main():
    args = parse_args()
    cyto = load_cyto(args.cyto)

    # chromosome sizes
    chrom_size = cyto.groupby("chrom")["end"].max().to_dict()

    # chromosomes to draw
    chroms = human_chrom_order(list(chrom_size.keys()))
    CHROMS = chroms

    # optional KZFP/genes
    kzfp = None
    if args.bed and os.path.isfile(args.bed):
        kzfp = pd.read_csv(args.bed, sep="\t", header=None, names=["chrom","start","end","name"])
        kzfp["mid"] = ((kzfp["start"] + kzfp["end"])//2).astype(int)
        # 只留图上有的染色体
        kzfp = kzfp[kzfp["chrom"].isin(CHROMS)].copy()

    # layout
    cols = args.cols
    rows = int(math.ceil(len(CHROMS)/cols))
    fig_w, fig_h = 16, 12
    W = 0.55         # chr width (x direction)
    XGAP = 2.2       # column gap
    YGAP = 1.0       # row gap
    TOP = rows*YGAP + 0.8  # top y

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")

    chrom_box = {}  # chrom -> (x0, y0)

    for i, chrom in enumerate(CHROMS):
        col = i % cols
        row = i // cols
        x0 = col*XGAP
        y0 = TOP - row*YGAP

        chrom_box[chrom] = (x0, y0)

        L = chrom_size[chrom] * 1.0  # length in bp

        # chromosome backbone (outline)
        ax.add_patch(patches.Rectangle((x0, y0), W, 1.0, fill=False, linewidth=0.8))

        # paint bands
        sub = cyto[cyto["chrom"]==chrom].sort_values("start")
        for _, r in sub.iterrows():
            s, e, stain = int(r["start"]), int(r["end"]), r["gieStain"]
            frac_s = s / L
            frac_e = e / L
            h = max(1e-6, frac_e - frac_s)
            face = band_color(stain)
            # draw band as filled rect
            ax.add_patch(patches.Rectangle((x0, y0+frac_s), W, h, facecolor=face, edgecolor="none"))

        # chromosome label
        ax.text(x0 + W + 0.15, y0 + 0.5, chrom.replace("chr",""), fontsize=9, va="center")

    # place labels if BED provided
    if kzfp is not None and len(kzfp):
        y_taken = defaultdict(list)
        for _, g in kzfp.sort_values(["chrom","mid","name"]).iterrows():
            chrom, pos, name = g["chrom"], int(g["mid"]), str(g["name"])
            if chrom not in chrom_box or chrom not in chrom_size:
                continue
            L = chrom_size[chrom]*1.0
            frac = np.clip(pos / L, 0, 1)
            x0, y0 = chrom_box[chrom]
            y = y0 + frac
            x_tick = x0 + W

            # avoid overlap (simple)
            dy = 0.0
            step = 0.02
            while any(abs((y+dy) - yy) < 0.02 for yy in y_taken[chrom]):
                dy += step if (len(y_taken[chrom]) % 2 == 0) else -step
            y_final = y + dy
            y_taken[chrom].append(y_final)

            # lead line
            ax.plot([x_tick, x_tick+0.6], [y, y_final], linewidth=0.6)
            # label text with band name
            band = label_for_band(cyto, chrom, pos)
            txt = f"{name} {band}" if band else name
            ax.text(x_tick+0.65, y_final, txt, fontsize=6, va="center")

    # set limits to ensure everything visible
    ax.set_xlim(-0.5, cols*XGAP - (XGAP - (W+1.4)))  # leave room for labels
    ax.set_ylim(-0.2, TOP + 0.2)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print("Saved:", args.out)

if __name__ == "__main__":
    main()
