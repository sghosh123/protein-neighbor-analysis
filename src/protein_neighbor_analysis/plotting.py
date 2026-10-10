"""Optional plotting helpers for frame summaries and residue occupancy."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def _plt():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Plotting requires matplotlib; install with pip install '.[plot]'") from exc


def plot_json_summaries(frame_rows: list[dict[str, Any]], occupancy_rows: list[dict[str, Any]], output_dir: str | Path) -> list[Path]:
    """Write per-source neighbor-count traces and residue-occupancy bars."""
    from collections import defaultdict

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    plt = _plt()
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in frame_rows:
        grouped[str(row["source_file"])].append(row)
    written: list[Path] = []
    for source, rows in grouped.items():
        rows = sorted(rows, key=lambda r: r["frame_ordinal_within_file"])
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot([r["frame_ordinal_within_file"] for r in rows], [r["n_unique_neighbor_atoms"] for r in rows], marker=".", linewidth=1)
        ax.set_xlabel("Frame ordinal within source file")
        ax.set_ylabel("Unique neighbor atoms")
        ax.set_title(source)
        fig.tight_layout()
        target = out / f"{Path(source).stem}_neighbor_trace.png"
        fig.savefig(target, dpi=160)
        plt.close(fig)
        written.append(target)

        counts = [int(r["n_unique_neighbor_atoms"]) for r in rows]
        min_count, max_count = min(counts), max(counts)
        edges = [n - 0.5 for n in range(min_count, max_count + 2)]
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ax.hist(counts, bins=edges, edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Unique neighbor atoms in frame")
        ax.set_ylabel("Number of frames")
        ax.set_title(f"Neighbor-count distribution: {source}")
        fig.tight_layout()
        target = out / f"{Path(source).stem}_neighbor_count_distribution.png"
        fig.savefig(target, dpi=160)
        plt.close(fig)
        written.append(target)

        source_occ = [r for r in occupancy_rows if r["source_file"] == source]
        if source_occ:
            source_occ.sort(key=lambda r: (-float(r["occupancy_fraction"]), str(r["segid"]), str(r["resid"]), str(r["resname"])))
            # Limit only display complexity; the full occupancy table is retained.
            visible = source_occ[:30]
            labels = [f'{r["segid"]}:{r["resname"]}{r["resid"]}' for r in visible]
            fig, ax = plt.subplots(figsize=(max(8, len(visible) * 0.28), 4.5))
            ax.bar(range(len(visible)), [float(r["occupancy_fraction"]) for r in visible])
            ax.set_xticks(range(len(visible)), labels, rotation=75, ha="right", fontsize=8)
            ax.set_ylabel("Frame occupancy fraction")
            ax.set_title(f"Residue occupancy: {source} (top 30 shown at most)")
            fig.tight_layout()
            target = out / f"{Path(source).stem}_residue_occupancy.png"
            fig.savefig(target, dpi=160)
            plt.close(fig)
            written.append(target)
    return written
