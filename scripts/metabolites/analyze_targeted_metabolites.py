#!/usr/bin/env python3
"""
Metabolomics analysis for Pseudomonas targeted metabolites.

Steps:
  1. Filter: max(experimental) > 3x max(control)
  2. PCoA on filtered data (inoc + d0-non samples)
  3. Volcano plots: inoc vs non (media control, same timepoint)
  4. Volcano plots: inoc vs d0 (pre-incubation baseline)
  5. Additional plots: heatmap, time-course boxplots
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────────
CSV_PATH = "/Users/mingfeichen/Recyc_necromass_metabolites/Pseudomonas_targeted_metabolites.csv"
OUT_DIR  = os.path.dirname(CSV_PATH)

# ── helper: parse one column name ──────────────────────────────────────────
def parse_col(col):
    m = re.match(r"^\d+_(.+)_(\d+)$", col)
    if not m:
        return None
    name, rep = m.group(1), m.group(2)

    if "ExCtrl" in name:
        return {"col": col, "type": "ctrl", "group": "ExCtrl", "replicate": rep}

    parts = name.split("-")

    # noC controls
    if "noC" in name:
        return {"col": col, "type": "ctrl", "group": "noC", "replicate": rep}

    # TxCtrl d0  (e.g. TxCtrl-100C-Pseu-A-fresh-d0-non)
    if parts[0] == "TxCtrl":
        return {
            "col": col, "type": "TxCtrl_d0",
            "carbon": parts[1], "organism": parts[2],
            "aer_ana": parts[3], "fresh_recyc": parts[4],
            "day": parts[5], "inoc_non": parts[6],
            "replicate": rep,
        }

    # sup- experimental samples
    if parts[0] == "sup":
        return {
            "col": col, "type": "sup",
            "carbon": parts[1], "organism": parts[2],
            "aer_ana": parts[3], "fresh_recyc": parts[4],
            "day": parts[5], "inoc_non": parts[6],
            "replicate": rep,
        }
    return None


# ── load ───────────────────────────────────────────────────────────────────
print("Loading data …")
df = pd.read_csv(CSV_PATH, index_col=0)
df = df.apply(pd.to_numeric, errors="coerce")
print(f"  {df.shape[0]} metabolites × {df.shape[1]} samples")

meta = {c: parse_col(c) for c in df.columns}
meta = {k: v for k, v in meta.items() if v is not None}

ctrl_cols = [c for c, v in meta.items() if v["type"] == "ctrl"]
exp_cols  = [c for c in df.columns if c in meta and meta[c]["type"] != "ctrl"]

print(f"  ctrl cols : {len(ctrl_cols)}")
print(f"  exp cols  : {len(exp_cols)}")


# ── step 1: filter rows ────────────────────────────────────────────────────
print("\nFiltering rows (max_exp > 3 × max_ctrl) …")
max_ctrl = df[ctrl_cols].max(axis=1, skipna=True).fillna(0)
max_exp  = df[exp_cols].max(axis=1, skipna=True)
mask = (max_exp > 3 * max_ctrl) & (max_exp > 0)
df_f = df.loc[mask].copy()
print(f"  {mask.sum()} / {len(mask)} metabolites pass")


# ── imputation ─────────────────────────────────────────────────────────────
def impute(df_in):
    out = df_in.copy().astype(float)
    for idx in out.index:
        row = out.loc[idx]
        pos = row[row > 0]
        fv = pos.min() / 2 if len(pos) else 1.0
        out.loc[idx] = row.fillna(fv)
    return out


# ── volcano helper ─────────────────────────────────────────────────────────
def make_volcano(df_imp, grp_a, grp_b, title, out_stem,
                 xlabel="log₂FC (inoc / reference)",
                 fc_thr=1.0, p_thr=0.05, label_n=6):
    """Welch t-test on log2-values; volcano with top-hit labels."""
    rows = []
    for met in df_imp.index:
        a = df_imp.loc[met, grp_a].dropna().values.astype(float)
        b = df_imp.loc[met, grp_b].dropna().values.astype(float)
        if len(a) == 0 or len(b) == 0:
            continue
        la, lb = np.log2(np.where(a == 0, 1e-6, a)), np.log2(np.where(b == 0, 1e-6, b))
        fc = la.mean() - lb.mean()
        pv = stats.ttest_ind(la, lb, equal_var=False).pvalue if (len(a) >= 2 and len(b) >= 2) else np.nan
        rows.append({"metabolite": met, "log2FC": fc, "pval": pv})

    if not rows:
        print(f"  [skip] {title}: no data")
        return pd.DataFrame()

    res = pd.DataFrame(rows)
    res["-log10p"] = -np.log10(res["pval"].clip(lower=1e-300))
    sig_up   = (res["pval"] < p_thr) & (res["log2FC"] >  fc_thr)
    sig_down = (res["pval"] < p_thr) & (res["log2FC"] < -fc_thr)

    fig, ax = plt.subplots(figsize=(7, 6))

    # non-significant
    ax.scatter(res.loc[~sig_up & ~sig_down, "log2FC"],
               res.loc[~sig_up & ~sig_down, "-log10p"],
               c="#AAAAAA", s=35, alpha=0.6, linewidths=0, label="NS")

    def _label_top(sub_df, col="log2FC", n=label_n):
        top = sub_df.nlargest(n, "-log10p") if col == "log2FC" else sub_df.nsmallest(n, "log2FC")
        for _, row in top.iterrows():
            short = str(row["metabolite"])
            # extract readable name between first _ and last _xxx
            m2 = re.match(r"^\d+_(.+?)(?:_positive|_negative).*$", short)
            label = m2.group(1).replace("_", " ") if m2 else short.split("_")[1] if "_" in short else short
            ax.annotate(label, (row["log2FC"], row["-log10p"]),
                        fontsize=5.5, ha="left" if col == "log2FC" else "right",
                        xytext=(3, 2), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color="black", lw=0.4))

    if sig_up.sum():
        ax.scatter(res.loc[sig_up, "log2FC"], res.loc[sig_up, "-log10p"],
                   c="#E74C3C", s=55, alpha=0.85, linewidths=0,
                   label=f"Up in inoc (n={sig_up.sum()})")
        _label_top(res[sig_up], "log2FC")

    if sig_down.sum():
        ax.scatter(res.loc[sig_down, "log2FC"], res.loc[sig_down, "-log10p"],
                   c="#2980B9", s=55, alpha=0.85, linewidths=0,
                   label=f"Down in inoc (n={sig_down.sum()})")
        _label_top(res[sig_down], "wrong_col")   # pass smallest log2FC

    ax.axhline(-np.log10(p_thr), color="k", lw=0.8, ls="--", alpha=0.5)
    ax.axvline( fc_thr, color="k", lw=0.8, ls="--", alpha=0.5)
    ax.axvline(-fc_thr, color="k", lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=10)
    ax.set_title(title, fontsize=9, pad=6)
    ax.legend(fontsize=7, loc="upper left")
    plt.tight_layout()
    for ext in ("pdf", "png"):
        plt.savefig(f"{out_stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close()
    return res


# ── select experimental columns ────────────────────────────────────────────
df_exp = df_f[exp_cols].copy()
df_imp = impute(df_exp)

sup_meta    = {c: meta[c] for c in exp_cols if meta[c]["type"] == "sup"}
txctrl_meta = {c: meta[c] for c in exp_cols if meta[c]["type"] == "TxCtrl_d0"}


# ══════════════════════════════════════════════════════════════════════════
# Step 2: PCoA (PCA on log2-standardised data)
# ══════════════════════════════════════════════════════════════════════════
print("\nRunning PCoA …")

inoc_cols = [c for c, v in sup_meta.items() if v["inoc_non"] == "inoc"]
d0_cols   = list(txctrl_meta.keys())
pcoa_cols = inoc_cols + d0_cols

X_pcoa = np.log2(df_imp[pcoa_cols] + 1).T.values   # (samples, features)
X_pcoa = StandardScaler().fit_transform(X_pcoa)

pca  = PCA(n_components=min(5, X_pcoa.shape[0] - 1))
pcs  = pca.fit_transform(X_pcoa)

# build per-sample metadata
pcoa_rows = []
for i, c in enumerate(pcoa_cols):
    if c in sup_meta:
        v = sup_meta[c]
        pcoa_rows.append({
            "PC1": pcs[i, 0], "PC2": pcs[i, 1],
            "aer_ana": v["aer_ana"],
            "fresh_recyc": v["fresh_recyc"],
            "day": v["day"],
            "inoc_non": "inoc",
        })
    else:
        v = txctrl_meta[c]
        pcoa_rows.append({
            "PC1": pcs[i, 0], "PC2": pcs[i, 1],
            "aer_ana": v["aer_ana"],
            "fresh_recyc": v["fresh_recyc"],
            "day": "d0",
            "inoc_non": "non",
        })

pm = pd.DataFrame(pcoa_rows)

# colours / shapes
col_map   = {"A": "#2196F3", "Ana": "#E64980"}
day_order = ["d0", "d2", "d6", "d11", "d27"]
mrk_map   = dict(zip(day_order, ["s", "o", "^", "D", "*"]))
edge_map  = {"fresh": "#222222", "recyc": "#AAAAAA"}

fig, ax = plt.subplots(figsize=(9, 7))

for _, row in pm.iterrows():
    c    = "#4CAF50" if row["inoc_non"] == "non" else col_map.get(row["aer_ana"], "#999999")
    mrk  = mrk_map.get(row["day"], "o")
    ec   = edge_map.get(row["fresh_recyc"], "black")
    sz   = 130 if row["inoc_non"] == "non" else 90
    ax.scatter(row["PC1"], row["PC2"], c=c, marker=mrk, edgecolors=ec,
               linewidths=1.4, s=sz, alpha=0.85, zorder=3)

# legend
leg = []
leg.append(mpatches.Patch(fc="#4CAF50",  label="d0 non-inoculated (baseline)"))
leg.append(mpatches.Patch(fc="#2196F3",  label="Aerobic inoculated"))
leg.append(mpatches.Patch(fc="#E64980",  label="Anaerobic inoculated"))
for d in day_order:
    leg.append(Line2D([0], [0], marker=mrk_map[d], color="w",
                      markerfacecolor="#777777", markeredgecolor="k",
                      markersize=8, label=d))
leg.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#777",
                  markeredgecolor="#222", markersize=8, label="fresh"))
leg.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#777",
                  markeredgecolor="#AAA", markersize=8, label="recycled"))

ax.legend(handles=leg, fontsize=7.5, loc="upper right",
          framealpha=0.9, ncol=2)
ax.axhline(0, color="gray", lw=0.5, ls="--", alpha=0.4)
ax.axvline(0, color="gray", lw=0.5, ls="--", alpha=0.4)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f} %)", fontsize=11)
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f} %)", fontsize=11)
ax.set_title("PCoA of filtered metabolites\n"
             "(colour = aer/ana, shape = day, edge = fresh/recycled)", fontsize=10)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(OUT_DIR, f"pcoa.{ext}"), dpi=200, bbox_inches="tight")
plt.close()
print("  saved pcoa.pdf / pcoa.png")


# ══════════════════════════════════════════════════════════════════════════
# Step 3: Volcano plots — inoc vs non (same timepoint, same condition)
# ══════════════════════════════════════════════════════════════════════════
print("\nVolcano plots (inoc vs non) …")
vdir3 = os.path.join(OUT_DIR, "volcano_inoc_vs_non")
os.makedirs(vdir3, exist_ok=True)

# collect unique (aer_ana, fresh_recyc, day) that have BOTH inoc AND non
conditions_inoc = set((v["aer_ana"], v["fresh_recyc"], v["day"])
                      for v in sup_meta.values() if v["inoc_non"] == "inoc")
conditions_non  = set((v["aer_ana"], v["fresh_recyc"], v["day"])
                      for v in sup_meta.values() if v["inoc_non"] == "non")
paired = sorted(conditions_inoc & conditions_non)

all_v3 = {}
for (aer, fr, day) in paired:
    g_inoc = [c for c, v in sup_meta.items()
              if v["aer_ana"]==aer and v["fresh_recyc"]==fr
              and v["day"]==day and v["inoc_non"]=="inoc"]
    g_non  = [c for c, v in sup_meta.items()
              if v["aer_ana"]==aer and v["fresh_recyc"]==fr
              and v["day"]==day and v["inoc_non"]=="non"]
    stem  = os.path.join(vdir3, f"volcano_{aer}_{fr}_{day}_inoc_vs_non")
    title = f"{aer} · {fr} · {day}  —  inoc vs non"
    res = make_volcano(df_imp, g_inoc, g_non, title, stem,
                       xlabel="log₂FC (inoc / non)")
    if not res.empty:
        res["condition"] = f"{aer}-{fr}-{day}"
        all_v3[f"{aer}_{fr}_{day}"] = res
    print(f"  {title}: {(res['pval']<0.05).sum() if not res.empty else 0} sig.")

# save combined table
if all_v3:
    pd.concat(all_v3.values()).to_csv(
        os.path.join(OUT_DIR, "inoc_vs_non_results.csv"), index=False)
    print("  saved inoc_vs_non_results.csv")


# ══════════════════════════════════════════════════════════════════════════
# Step 4: Volcano plots — inoc vs d0 baseline
# ══════════════════════════════════════════════════════════════════════════
print("\nVolcano plots (inoc vs d0) …")
vdir4 = os.path.join(OUT_DIR, "volcano_inoc_vs_d0")
os.makedirs(vdir4, exist_ok=True)

all_v4 = {}
for (aer, fr, day) in sorted(conditions_inoc):
    g_inoc = [c for c, v in sup_meta.items()
              if v["aer_ana"]==aer and v["fresh_recyc"]==fr
              and v["day"]==day and v["inoc_non"]=="inoc"]
    # d0 cols matching fresh_recyc
    g_d0 = [c for c, v in txctrl_meta.items()
             if v["fresh_recyc"] == fr]
    if not g_d0:
        g_d0 = list(txctrl_meta.keys())   # fallback: all d0
    stem  = os.path.join(vdir4, f"volcano_{aer}_{fr}_{day}_inoc_vs_d0")
    title = f"{aer} · {fr} · {day}  —  inoc vs d0"
    res = make_volcano(df_imp, g_inoc, g_d0, title, stem,
                       xlabel="log₂FC (inoc / d0)")
    if not res.empty:
        res["condition"] = f"{aer}-{fr}-{day}"
        all_v4[f"{aer}_{fr}_{day}"] = res
    print(f"  {title}: {(res['pval']<0.05).sum() if not res.empty else 0} sig.")

if all_v4:
    pd.concat(all_v4.values()).to_csv(
        os.path.join(OUT_DIR, "inoc_vs_d0_results.csv"), index=False)
    print("  saved inoc_vs_d0_results.csv")


# ══════════════════════════════════════════════════════════════════════════
# Step 5a: Heatmap of filtered metabolites (aerobic inoc + d0 samples)
# ══════════════════════════════════════════════════════════════════════════
print("\nHeatmap …")

hmap_cols = [c for c, v in sup_meta.items()
             if v["aer_ana"] == "A" and v["inoc_non"] == "inoc"] + d0_cols
hmap_cols = [c for c in hmap_cols if c in df_imp.columns]

# column order: d0 first, then by day
def sort_key(c):
    if c in txctrl_meta:
        v = txctrl_meta[c]
        return (0, v["fresh_recyc"], v["replicate"])
    v = sup_meta[c]
    day_n = int(v["day"][1:]) if v["day"][1:].isdigit() else 0
    return (day_n, v["fresh_recyc"], v["inoc_non"], v["replicate"])

hmap_cols_sorted = sorted(hmap_cols, key=sort_key)

X_h  = np.log2(df_imp[hmap_cols_sorted] + 1)
row_z = X_h.subtract(X_h.mean(axis=1), axis=0).divide(X_h.std(axis=1).replace(0, 1), axis=0)

# short metabolite labels
def short_name(s):
    m2 = re.match(r"^\d+_(.+?)(?:_positive|_negative).*$", s)
    if m2:
        return m2.group(1).replace("_", " ")[:30]
    return s[:30]

row_labels = [short_name(r) for r in row_z.index]

col_labels = []
for c in hmap_cols_sorted:
    if c in txctrl_meta:
        v = txctrl_meta[c]
        col_labels.append(f"d0-{v['fresh_recyc'][:2]}")
    else:
        v = sup_meta[c]
        col_labels.append(f"{v['day']}-{v['fresh_recyc'][:2]}")

# cluster rows
from scipy.cluster import hierarchy
link = hierarchy.linkage(row_z.fillna(0), method="ward", metric="euclidean")
order = hierarchy.dendrogram(link, no_plot=True)["leaves"]
row_z_cl = row_z.iloc[order]
row_labels_cl = [row_labels[i] for i in order]

fig_h, ax_h = plt.subplots(figsize=(max(12, len(hmap_cols_sorted)*0.35),
                                     max(8, len(row_z)*0.18)))
im = ax_h.imshow(row_z_cl.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax_h.set_xticks(range(len(hmap_cols_sorted)))
ax_h.set_xticklabels(col_labels, rotation=90, fontsize=6)
ax_h.set_yticks(range(len(row_z_cl)))
ax_h.set_yticklabels(row_labels_cl, fontsize=6)
plt.colorbar(im, ax=ax_h, label="z-score (log₂ intensity)", shrink=0.6)
ax_h.set_title("Heatmap of filtered metabolites (aerobic inoc + d0)", fontsize=10)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(OUT_DIR, f"heatmap.{ext}"), dpi=200, bbox_inches="tight")
plt.close()
print("  saved heatmap.pdf / heatmap.png")


# ══════════════════════════════════════════════════════════════════════════
# Step 5b: Time-course box plots for top enriched metabolites
#          (pick top metabolites by overall variance across inoc samples)
# ══════════════════════════════════════════════════════════════════════════
print("\nTime-course box plots …")

inoc_cols_A = [c for c, v in sup_meta.items()
               if v["aer_ana"] == "A" and v["inoc_non"] == "inoc"]

day_order_num = ["d0", "d2", "d6", "d11", "d27"]

# variance across aerobic inoc columns
var_scores = df_imp[inoc_cols_A].var(axis=1)
top_mets   = var_scores.nlargest(12).index.tolist()

fig_tc, axes = plt.subplots(3, 4, figsize=(18, 12))
axes = axes.flatten()

for ax, met in zip(axes, top_mets):
    # collect data by day and fresh/recyc
    plot_data = {}
    for fr in ["fresh", "recyc"]:
        for day in day_order_num:
            if day == "d0":
                cols = [c for c, v in txctrl_meta.items() if v["fresh_recyc"] == fr]
            else:
                cols = [c for c, v in sup_meta.items()
                        if v["aer_ana"] == "A" and v["fresh_recyc"] == fr
                        and v["day"] == day and v["inoc_non"] == "inoc"]
            if cols:
                vals = np.log2(df_imp.loc[met, cols].values.astype(float) + 1)
                plot_data[(fr, day)] = vals

    # plot
    x_pos, x_ticks, x_labels = 0, [], []
    for fr, color in [("fresh", "#2196F3"), ("recyc", "#FF9800")]:
        for day in day_order_num:
            if (fr, day) in plot_data:
                vals = plot_data[(fr, day)]
                bp = ax.boxplot(vals, positions=[x_pos], widths=0.6,
                                patch_artist=True,
                                boxprops=dict(facecolor=color, alpha=0.6),
                                medianprops=dict(color="black", lw=1.5),
                                whiskerprops=dict(color=color),
                                capprops=dict(color=color),
                                flierprops=dict(marker=".", color=color, ms=4))
                x_ticks.append(x_pos)
                x_labels.append(day)
                x_pos += 1
        x_pos += 0.5   # gap between fresh and recyc

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=6, rotation=45)
    ax.set_ylabel("log₂ intensity", fontsize=7)
    short = short_name(met)
    ax.set_title(short, fontsize=7, pad=3)
    ax.tick_params(labelsize=6)

# shared legend
legend_patches = [
    mpatches.Patch(facecolor="#2196F3", alpha=0.6, label="fresh"),
    mpatches.Patch(facecolor="#FF9800", alpha=0.6, label="recycled"),
]
fig_tc.legend(handles=legend_patches, loc="lower right", fontsize=9)
fig_tc.suptitle("Time-course (aerobic inoculated) — top variable metabolites",
                fontsize=11, y=1.01)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(OUT_DIR, f"timecourse_top12.{ext}"),
                dpi=200, bbox_inches="tight")
plt.close()
print("  saved timecourse_top12.pdf / timecourse_top12.png")


# ══════════════════════════════════════════════════════════════════════════
# Step 5c: Summary bar chart — number of significant metabolites per contrast
# ══════════════════════════════════════════════════════════════════════════
print("\nSummary bar chart …")

def sig_counts(d, fc_thr=1, p_thr=0.05):
    out = {}
    for k, res in d.items():
        up   = ((res["log2FC"] >  fc_thr) & (res["pval"] < p_thr)).sum()
        down = ((res["log2FC"] < -fc_thr) & (res["pval"] < p_thr)).sum()
        out[k] = (up, down)
    return out

sc3 = sig_counts(all_v3)
sc4 = sig_counts(all_v4)

def _bar_panel(ax, sc, title):
    if not sc:
        ax.set_visible(False)
        return
    keys  = list(sc.keys())
    ups   = [sc[k][0] for k in keys]
    downs = [-sc[k][1] for k in keys]
    x     = np.arange(len(keys))
    ax.bar(x, ups,   color="#E74C3C", label="Up in inoc")
    ax.bar(x, downs, color="#2980B9", label="Down in inoc")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([k.replace("_", "\n") for k in keys], fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("# significant metabolites")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)

fig_s, (axA, axB) = plt.subplots(1, 2, figsize=(16, 5))
_bar_panel(axA, sc3, "Inoc vs Non (same timepoint)")
_bar_panel(axB, sc4, "Inoc vs d0 baseline")
plt.suptitle("Significant metabolites (|log₂FC|>1, p<0.05)", fontsize=11)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(OUT_DIR, f"summary_sig_counts.{ext}"),
                dpi=200, bbox_inches="tight")
plt.close()
print("  saved summary_sig_counts.pdf / summary_sig_counts.png")


# ══════════════════════════════════════════════════════════════════════════
# Step 5d: Aer vs Ana PCoA (separate panel)
# ══════════════════════════════════════════════════════════════════════════
print("\nAer vs Ana PCoA …")
ana_inoc = [c for c, v in sup_meta.items() if v["inoc_non"] == "inoc"]
X2 = np.log2(df_imp[ana_inoc] + 1).T.values
X2 = StandardScaler().fit_transform(X2)
pca2 = PCA(n_components=2)
pcs2 = pca2.fit_transform(X2)

rows2 = []
for i, c in enumerate(ana_inoc):
    v = sup_meta[c]
    rows2.append({"PC1": pcs2[i,0], "PC2": pcs2[i,1],
                  "aer_ana": v["aer_ana"], "fresh_recyc": v["fresh_recyc"],
                  "day": v["day"]})
pm2 = pd.DataFrame(rows2)

fig2, ax2 = plt.subplots(figsize=(8, 6))
for _, row in pm2.iterrows():
    c    = "#2196F3" if row["aer_ana"] == "A" else "#E64980"
    mrk  = mrk_map.get(row["day"], "o")
    ec   = edge_map.get(row["fresh_recyc"], "black")
    ax2.scatter(row["PC1"], row["PC2"], c=c, marker=mrk, edgecolors=ec,
                linewidths=1.3, s=90, alpha=0.85, zorder=3)

leg2 = [mpatches.Patch(fc="#2196F3", label="Aerobic"),
        mpatches.Patch(fc="#E64980", label="Anaerobic")]
for d in ["d2", "d6", "d11", "d27"]:
    leg2.append(Line2D([0], [0], marker=mrk_map[d], color="w",
                       markerfacecolor="#777", markeredgecolor="k",
                       markersize=8, label=d))
leg2.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#777",
                   markeredgecolor="#222", markersize=8, label="fresh"))
leg2.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#777",
                   markeredgecolor="#AAA", markersize=8, label="recycled"))
ax2.legend(handles=leg2, fontsize=7.5, loc="upper right", ncol=2)
ax2.axhline(0, color="gray", lw=0.5, ls="--", alpha=0.4)
ax2.axvline(0, color="gray", lw=0.5, ls="--", alpha=0.4)
ax2.set_xlabel(f"PC1 ({pca2.explained_variance_ratio_[0]*100:.1f} %)", fontsize=11)
ax2.set_ylabel(f"PC2 ({pca2.explained_variance_ratio_[1]*100:.1f} %)", fontsize=11)
ax2.set_title("PCoA — inoculated samples only (aer vs ana)", fontsize=10)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(os.path.join(OUT_DIR, f"pcoa_aer_vs_ana.{ext}"),
                dpi=200, bbox_inches="tight")
plt.close()
print("  saved pcoa_aer_vs_ana.pdf / pcoa_aer_vs_ana.png")


print("\n✓ All done. Output folder:", OUT_DIR)
