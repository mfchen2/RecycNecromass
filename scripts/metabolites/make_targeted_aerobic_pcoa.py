#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform


ROOT = Path("/Users/mingfeichen/Recyc_necromass")
OUTSETS = {
    "pseudomonas": {
        "csv": ROOT / "outputs" / "targeted_metabolites_filtered" / "Pseudo_targeted_metabolites_051926_filtered_3x.csv",
        "outdir": ROOT / "outputs" / "targeted_metabolites_pseudomonas",
        "title": "Pseudomonas targeted metabolites PCoA",
    },
    "arthrobacter": {
        "csv": ROOT / "outputs" / "targeted_metabolites_filtered" / "Arthrobacter_targeted_metabolites_filtered_3x.csv",
        "outdir": ROOT / "outputs" / "targeted_metabolites_arthrobacter",
        "title": "Arthrobacter targeted metabolites PCoA",
    },
}

DAY_ORDER_WITH_D0 = ["d0", "d2", "d6", "d11", "d27"]
DAY_ORDER_NO_D0 = ["d2", "d6", "d11", "d27"]
DAY_MARKERS = {"d0": "s", "d2": "o", "d6": "^", "d11": "D", "d27": "*"}
BLUE = "#1F77B4"
GREEN = "#1B9E77"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "legend.title_fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
        }
    )


def parse_sample_column(col: str) -> dict | None:
    m = re.match(r"^\d+_(.+)_(\d+)$", col)
    if not m:
        return None
    name, rep = m.group(1), m.group(2)
    parts = name.split("-")
    if parts[0] == "TxCtrl":
        if len(parts) < 7:
            return None
        return {
            "column": col,
            "sample_id": name,
            "replicate": int(rep),
            "type": "d0",
            "organism": parts[2],
            "oxygen": parts[3],
            "source": parts[4],
            "day": parts[5],
            "inoculation": parts[6],
        }
    if parts[0] == "sup":
        if len(parts) < 7:
            return None
        return {
            "column": col,
            "sample_id": name,
            "replicate": int(rep),
            "type": "sup",
            "organism": parts[2],
            "oxygen": parts[3],
            "source": parts[4],
            "day": parts[5],
            "inoculation": parts[6],
        }
    return None


def pcoa_from_bray(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dist = squareform(pdist(matrix, metric="braycurtis"))
    n = dist.shape[0]
    d2 = dist ** 2
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ d2 @ h
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 0
    eigvals_pos = eigvals[positive]
    eigvecs_pos = eigvecs[:, positive]
    if len(eigvals_pos) < 2:
        raise ValueError("Not enough positive eigenvalues for 2D PCoA.")
    coords = eigvecs_pos[:, :2] * np.sqrt(eigvals_pos[:2])
    pct = 100 * eigvals_pos[:2] / eigvals_pos.sum()
    return coords, np.round(pct, 1)


def load_panel(csv_path: Path, aerobic_only: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(csv_path)
    sample_meta = []
    sample_cols = []
    for col in df.columns:
        parsed = parse_sample_column(col)
        if parsed is not None:
            sample_meta.append(parsed)
            sample_cols.append(col)
    meta = pd.DataFrame(sample_meta)
    if meta.empty:
        raise ValueError(f"No sample columns found in {csv_path}")

    if aerobic_only:
        meta = meta.loc[meta["oxygen"].astype(str).eq("A")].copy()
    meta = meta.sort_values(["type", "source", "day", "replicate"]).reset_index(drop=True)
    meta["source"] = pd.Categorical(meta["source"], categories=["fresh", "recyc"], ordered=True)
    meta["day"] = pd.Categorical(meta["day"], categories=DAY_ORDER_WITH_D0, ordered=True)

    cols = meta["column"].tolist()
    mat = df.loc[:, cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).T.to_numpy()
    mat = np.log1p(mat)
    keep_features = mat.sum(axis=0) > 0
    mat = mat[:, keep_features]
    if mat.shape[1] < 2:
        raise ValueError(f"Not enough nonzero features in {csv_path}")

    coords, pct = pcoa_from_bray(mat)
    scores = meta.copy()
    scores["PCoA1"] = coords[:, 0]
    scores["PCoA2"] = coords[:, 1]
    scores.attrs["pct"] = pct
    return scores, df


def split_scores(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pct = scores.attrs.get("pct")
    with_d0_parts = [
        scores.loc[scores["type"] == "d0"].copy().reset_index(drop=True),
        scores.loc[(scores["type"] == "sup") & (scores["inoculation"].astype(str) == "inoc")].copy().reset_index(drop=True),
    ]
    for part in with_d0_parts:
        part.attrs = {}
    with_d0 = pd.concat(with_d0_parts, ignore_index=True)
    with_d0.attrs["pct"] = pct
    no_d0 = scores.loc[
        (scores["type"] == "sup")
        & (scores["inoculation"].astype(str) == "inoc")
        & (scores["day"].astype(str) != "d0")
    ].copy()
    no_d0.attrs = {}
    with_d0 = with_d0.sort_values(["type", "source", "day", "replicate"]).reset_index(drop=True)
    no_d0 = no_d0.sort_values(["source", "day", "replicate"]).reset_index(drop=True)
    no_d0.attrs["pct"] = pct
    return with_d0, no_d0


def plot_pcoa(scores: pd.DataFrame, title: str, out_base: Path, include_d0: bool) -> None:
    pct = scores.attrs["pct"]
    day_order = DAY_ORDER_WITH_D0 if include_d0 else DAY_ORDER_NO_D0

    fig, ax = plt.subplots(figsize=(8.7, 6.3))

    # experimental aerobic samples
    if include_d0:
        exp = scores.loc[(scores["type"] == "sup") & (scores["inoculation"].astype(str) == "inoc")].copy()
    else:
        exp = scores.loc[
            (scores["type"] == "sup")
            & (scores["inoculation"].astype(str) == "inoc")
            & (scores["day"].astype(str) != "d0")
        ].copy()

    for day in day_order:
        for source, fill in [("fresh", "none"), ("recyc", BLUE)]:
            sub = exp.loc[(exp["day"].astype(str) == day) & (exp["source"].astype(str) == source)]
            if sub.empty:
                continue
            ax.scatter(
                sub["PCoA1"],
                sub["PCoA2"],
                s=120,
                marker=DAY_MARKERS[day],
                facecolors=fill,
                edgecolors=BLUE,
                linewidths=1.2,
                alpha=0.95,
                zorder=3,
            )

    # d0 baseline samples
    if include_d0:
        d0 = scores.loc[scores["type"] == "d0"].copy()
        for source, fill in [("fresh", "none"), ("recyc", GREEN)]:
            sub = d0.loc[d0["source"].astype(str) == source]
            if sub.empty:
                continue
            ax.scatter(
                sub["PCoA1"],
                sub["PCoA2"],
                s=130,
                marker="s",
                facecolors=fill,
                edgecolors=GREEN,
                linewidths=1.2,
                alpha=0.95,
                zorder=4,
            )

    ax.axhline(0, color="#cfcfcf", lw=0.9, ls="--", zorder=0)
    ax.axvline(0, color="#cfcfcf", lw=0.9, ls="--", zorder=0)
    ax.set_xlabel(f"PCoA1 ({pct[0]}% of variance)")
    ax.set_ylabel(f"PCoA2 ({pct[1]}% of variance)")
    ax.set_title(title, fontsize=10.5, pad=8)
    ax.tick_params(direction="out", length=4, width=0.8)

    # legend blocks
    color_handles = [
        Line2D([0], [0], color=BLUE, lw=8, label="Aerobic"),
    ]
    if include_d0:
        color_handles.append(Line2D([0], [0], color=GREEN, lw=8, label="d0 baseline"))
    day_handles = [
        Line2D([0], [0], marker=DAY_MARKERS[day], color="black", linestyle="None", markersize=9, markerfacecolor="none" if day != "d27" else "black", markeredgewidth=1.0, label=day)
        for day in day_order
    ]
    source_handles = [
        Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=8, markerfacecolor="none", markeredgewidth=1.0, label="fresh (open)"),
        Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=8, markerfacecolor="black", markeredgewidth=1.0, label="recycled (filled)"),
    ]

    leg1 = ax.legend(handles=color_handles, title="Color", frameon=False, loc="upper right", bbox_to_anchor=(1.02, 1.00))
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=day_handles, title="Day", frameon=False, loc="center right", bbox_to_anchor=(1.02, 0.55))
    ax.add_artist(leg2)
    ax.legend(handles=source_handles, title="Source fill", frameon=False, loc="lower right", bbox_to_anchor=(1.02, 0.02))

    fig.tight_layout(rect=[0, 0, 0.86, 1])
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def process_one(name: str, spec: dict[str, Path | str]) -> None:
    scores, _ = load_panel(Path(spec["csv"]), aerobic_only=True)
    outdir = Path(spec["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    with_d0, no_d0 = split_scores(scores)

    with_d0.to_csv(outdir / "pcoa_targeted_metabolites_scores.csv", index=False)
    no_d0.to_csv(outdir / "pcoa_targeted_metabolites_no_d0_scores.csv", index=False)

    plot_pcoa(with_d0, str(spec["title"]), outdir / "pcoa_targeted_metabolites", include_d0=True)
    plot_pcoa(no_d0, str(spec["title"]), outdir / "pcoa_targeted_metabolites_no_d0", include_d0=False)


def main() -> None:
    configure_style()
    for name, spec in OUTSETS.items():
        process_one(name, spec)


if __name__ == "__main__":
    main()
