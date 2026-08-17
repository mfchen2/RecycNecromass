#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


ROOT = Path("/Users/mingfeichen/Recyc_necromass")
OUTSETS = {
    "pseudomonas": ROOT / "outputs" / "untargeted_metabolites_pseudomonas",
    "arthrobacter": ROOT / "outputs" / "untargeted_metabolites_arthrobacter",
}

SOURCE_COLORS = {"fresh": "#F8766D", "recyc": "#00BFC4"}
DAY_SHAPES = {"d0": "o", "d2": "s", "d6": "D", "d11": "^", "d27": "v"}
DAY_ORDER = ["d0", "d2", "d6", "d11", "d27"]
SOURCE_ORDER = ["fresh", "recyc"]


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.title_fontsize": 9,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_tables(out_dir: Path, sheet: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    values = pd.read_csv(out_dir / f"{sheet}_filtered_peak_heights.csv")
    meta = pd.read_csv(out_dir / f"{sheet}_sample_metadata.csv")
    return values, meta


def selected_aerobic_meta(meta: pd.DataFrame) -> pd.DataFrame:
    keep = meta["is_experimental"].astype(bool) & (
        ((meta["inoculation"] == "inoc") & (meta["day"].astype(str) != "d0"))
        | ((meta["inoculation"] == "non") & (meta["day"].astype(str) == "d0"))
    )
    keep &= meta["oxygen"].astype(str).eq("Aerobic")
    selected = meta.loc[keep].copy()
    selected["day"] = pd.Categorical(selected["day"].astype(str), categories=DAY_ORDER, ordered=True)
    selected["source"] = pd.Categorical(selected["source"].astype(str), categories=SOURCE_ORDER, ordered=True)
    return selected.sort_values(["source", "day", "replicate", "sample_id"]).reset_index(drop=True)


def classical_pcoa(distance_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = distance_matrix.shape[0]
    d2 = distance_matrix**2
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ d2 @ h
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 0
    eigvals_pos = eigvals[positive]
    eigvecs_pos = eigvecs[:, positive]
    coords = eigvecs_pos[:, :2] * np.sqrt(eigvals_pos[:2])
    return coords, eigvals_pos


def compute_pcoa(values: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    sample_cols = meta["column"].tolist()
    mat = values.loc[:, sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).T.to_numpy()
    mat = np.log1p(mat)
    mat = mat[:, mat.sum(axis=0) > 0]
    if mat.shape[1] < 2:
        raise ValueError("Not enough features after filtering for PCoA.")
    dist = squareform(pdist(mat, metric="braycurtis"))
    coords, eigvals_pos = classical_pcoa(dist)
    if coords.shape[1] < 2:
        raise ValueError("PCoA returned fewer than two axes.")
    total = eigvals_pos.sum()
    pct = np.round(100 * eigvals_pos[:2] / total, 1) if total > 0 else np.array([np.nan, np.nan])
    scores = meta.copy()
    scores["PCoA1"] = coords[:, 0]
    scores["PCoA2"] = coords[:, 1]
    scores.attrs["pct"] = pct
    return scores


def draw_pcoa(scores: pd.DataFrame, title: str, out_base: Path) -> None:
    pct = scores.attrs.get("pct", np.array([np.nan, np.nan]))
    fig, ax = plt.subplots(figsize=(8.6, 6.2))

    for day in DAY_ORDER:
        for source in SOURCE_ORDER:
            sub = scores[(scores["day"].astype(str) == day) & (scores["source"].astype(str) == source)]
            if sub.empty:
                continue
            ax.scatter(
                sub["PCoA1"],
                sub["PCoA2"],
                s=110,
                marker=DAY_SHAPES[day],
                c=SOURCE_COLORS[source],
                edgecolors=SOURCE_COLORS[source],
                linewidths=1.1,
                alpha=0.95,
                label=f"{source} | {day}",
            )

    ax.set_title(title, loc="left", fontweight="bold", pad=14)
    ax.set_xlabel(f"PCoA1 ({pct[0]}%)")
    ax.set_ylabel(f"PCoA2 ({pct[1]}%)")
    ax.axhline(0, color="#d9d9d9", lw=1, zorder=0)
    ax.axvline(0, color="#d9d9d9", lw=1, zorder=0)
    ax.grid(True, color="#ececec", linewidth=0.8)
    ax.set_axisbelow(True)

    from matplotlib.lines import Line2D

    day_handles = [
        Line2D(
            [0],
            [0],
            marker=DAY_SHAPES[day],
            color="black",
            linestyle="None",
            markersize=9,
            markerfacecolor="white",
            markeredgewidth=1.1,
            label=day,
        )
        for day in DAY_ORDER
    ]
    source_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            linestyle="None",
            markersize=10,
            markerfacecolor=SOURCE_COLORS[source],
            markeredgecolor=SOURCE_COLORS[source],
            label=source,
        )
        for source in SOURCE_ORDER
    ]

    legend1 = ax.legend(
        handles=day_handles,
        title="Day",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.66),
        borderaxespad=0.0,
        labelspacing=1.0,
        handletextpad=0.8,
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=source_handles,
        title="Source",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.26),
        borderaxespad=0.0,
        labelspacing=1.0,
        handletextpad=0.8,
    )

    fig.subplots_adjust(right=0.78, left=0.10, top=0.92, bottom=0.10)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def process_dataset(out_dir: Path) -> None:
    for sheet in ["positive", "negative"]:
        values, meta = read_tables(out_dir, sheet)
        selected = selected_aerobic_meta(meta)
        if selected.empty:
            continue
        scores = compute_pcoa(values, selected)
        scores.to_csv(out_dir / f"{sheet}_bray_pcoa_scores.csv", index=False)
        draw_pcoa(
            scores,
            f"{sheet} aerobic filtered features: Bray-Curtis PCoA",
            out_dir / f"{sheet}_bray_pcoa_selected_samples",
        )


def main() -> None:
    configure_matplotlib()
    for out_dir in OUTSETS.values():
        process_dataset(out_dir)


if __name__ == "__main__":
    main()
