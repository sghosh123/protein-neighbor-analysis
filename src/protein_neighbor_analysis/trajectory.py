"""Periodic-boundary-aware atom and residue neighborhood analysis."""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
except ImportError:  # pragma: no cover - availability checked at API boundary
    mda = None
    capped_distance = None

SUMMARY_FIELDS = [
    "frame", "time", "time_unit", "n_center_atoms", "n_contact_pairs",
    "n_neighbor_atoms", "n_neighbor_residues", "pbc_used",
    "box_x_A", "box_y_A", "box_z_A", "box_alpha_deg", "box_beta_deg", "box_gamma_deg",
]
CONTACT_FIELDS = [
    "frame", "time", "center_atom_index", "center_segid", "center_chain_id",
    "center_resid", "center_icode", "center_resname", "center_atom_name", "center_element", "center_resindex",
    "neighbor_atom_index", "neighbor_segid", "neighbor_chain_id", "neighbor_resid",
    "neighbor_icode", "neighbor_resname", "neighbor_atom_name", "neighbor_element",
    "distance_A", "neighbor_resindex",
]


def _safe_attr(obj: Any, attr: str, default: Any = "") -> Any:
    try:
        value = getattr(obj, attr)
        if value is None:
            return default
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            return default
        return value
    except (AttributeError, ValueError, TypeError):
        return default


def atom_record(atom: Any) -> dict[str, Any]:
    """Return stable, human-readable atom identifiers from an MDAnalysis atom."""
    return {
        "index": int(atom.index),
        "segid": str(_safe_attr(atom, "segid", "")),
        "chain_id": str(_safe_attr(atom, "chainID", "")),
        "resid": int(_safe_attr(atom, "resid", -1)),
        "icode": str(_safe_attr(atom, "icode", "")),
        "resname": str(_safe_attr(atom, "resname", "")),
        "name": str(_safe_attr(atom, "name", "")),
        "element": str(_safe_attr(atom, "element", "")),
        "resindex": int(_safe_attr(atom, "resindex", -1)),
    }


def valid_dimensions(dimensions: Any) -> bool:
    """True if dimensions are usable as an MDAnalysis unit cell."""
    if dimensions is None:
        return False
    try:
        arr = np.asarray(dimensions, dtype=float)
    except (TypeError, ValueError):
        return False
    if arr.shape != (6,) or not np.all(np.isfinite(arr)):
        return False
    return bool(np.all(arr[:3] > 0) and np.all((arr[3:] > 0) & (arr[3:] < 180)))


def _residue_identity(atom: Any) -> tuple[Any, ...]:
    """Use topology residue index to avoid merging residues with duplicate labels."""
    return (
        int(_safe_attr(atom, "resindex", -1)),
        str(_safe_attr(atom, "segid", "")),
        str(_safe_attr(atom, "chainID", "")),
        int(_safe_attr(atom, "resid", -1)),
        str(_safe_attr(atom, "icode", "")),
        str(_safe_attr(atom, "resname", "")),
    )


def analyze_universe(
    universe: Any,
    center_selection: str,
    search_selection: str = "all",
    cutoff: float = 3.5,
    *,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    pbc: str = "auto",
    exclude_center_atoms: bool = True,
    time_unit: str = "ps",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Analyze neighbors for each selected trajectory frame.

    Coordinates/distances are in Angstrom, as in MDAnalysis. PBC modes:
    ``auto`` uses the current frame box if valid; ``on`` requires a valid box;
    ``off`` calculates direct Cartesian distances. Contact rows represent atom
    pairs; frame summaries separately count pairs, unique neighboring atoms,
    and unique neighboring topology residues.
    """
    if capped_distance is None:
        raise ImportError("MDAnalysis is required. Install with: pip install '.[trajectory]'")
    if not isinstance(center_selection, str) or not center_selection.strip():
        raise ValueError("center_selection must be a non-empty MDAnalysis selection")
    if not isinstance(search_selection, str) or not search_selection.strip():
        raise ValueError("search_selection must be a non-empty MDAnalysis selection")
    if not math.isfinite(float(cutoff)) or cutoff <= 0:
        raise ValueError("cutoff must be a finite positive distance in Angstrom")
    if pbc not in ("auto", "on", "off"):
        raise ValueError("pbc must be one of: auto, on, off")
    if start < 0 or step <= 0:
        raise ValueError("start must be >= 0 and step must be > 0")
    n_frames = int(universe.trajectory.n_frames)
    if start >= n_frames:
        raise ValueError(f"start frame {start} is outside a {n_frames}-frame trajectory")
    effective_stop = n_frames if stop is None else min(int(stop), n_frames)
    if effective_stop <= start:
        raise ValueError("stop must be greater than start")

    # Parse selections exactly once; the AtomGroups update their positions by frame.
    try:
        center = universe.select_atoms(center_selection)
        search = universe.select_atoms(search_selection)
    except Exception as exc:
        raise ValueError(f"Invalid MDAnalysis selection: {exc}") from exc
    if len(center) == 0:
        raise ValueError(f"Center selection matched no atoms: {center_selection!r}")
    if len(search) == 0:
        raise ValueError(f"Search selection matched no atoms: {search_selection!r}")

    center_indices = set(map(int, center.indices))
    if exclude_center_atoms:
        keep = np.array([int(idx) not in center_indices for idx in search.indices], dtype=bool)
        search = search[keep]
        if len(search) == 0:
            raise ValueError("Search selection is empty after excluding center atoms")

    summary_rows: list[dict[str, Any]] = []
    contact_rows: list[dict[str, Any]] = []
    for ts in universe.trajectory[start:effective_stop:step]:
        dims = ts.dimensions
        box_ok = valid_dimensions(dims)
        if pbc == "on" and not box_ok:
            raise ValueError(f"PBC was required but frame {ts.frame} has no valid unit cell")
        box = np.asarray(dims, dtype=np.float32) if (pbc != "off" and box_ok) else None
        pairs, distances = capped_distance(
            center.positions,
            search.positions,
            max_cutoff=float(cutoff),
            box=box,
            return_distances=True,
        )
        frame = int(ts.frame)
        time_value = float(ts.time) if np.isfinite(ts.time) else float("nan")
        neighbor_local_indices: set[int] = set()
        neighbor_residue_ids: set[tuple[Any, ...]] = set()
        for pair, distance in zip(pairs, distances):
            ci, ni = int(pair[0]), int(pair[1])
            ca, na = center[ci], search[ni]
            c_rec, n_rec = atom_record(ca), atom_record(na)
            neighbor_local_indices.add(int(na.index))
            neighbor_residue_ids.add(_residue_identity(na))
            contact_rows.append({
                "frame": frame,
                "time": time_value,
                "center_atom_index": c_rec["index"],
                "center_segid": c_rec["segid"],
                "center_chain_id": c_rec["chain_id"],
                "center_resid": c_rec["resid"],
                "center_icode": c_rec["icode"],
                "center_resname": c_rec["resname"],
                "center_atom_name": c_rec["name"],
                "center_element": c_rec["element"],
                "center_resindex": c_rec["resindex"],
                "neighbor_atom_index": n_rec["index"],
                "neighbor_segid": n_rec["segid"],
                "neighbor_chain_id": n_rec["chain_id"],
                "neighbor_resid": n_rec["resid"],
                "neighbor_icode": n_rec["icode"],
                "neighbor_resname": n_rec["resname"],
                "neighbor_atom_name": n_rec["name"],
                "neighbor_element": n_rec["element"],
                "distance_A": float(distance),
                "neighbor_resindex": n_rec["resindex"],
            })
        summary_rows.append({
            "frame": frame,
            "time": time_value,
            "time_unit": time_unit,
            "n_center_atoms": len(center),
            "n_contact_pairs": len(pairs),
            "n_neighbor_atoms": len(neighbor_local_indices),
            "n_neighbor_residues": len(neighbor_residue_ids),
            "pbc_used": bool(box is not None),
            "box_x_A": float(dims[0]) if box_ok else "",
            "box_y_A": float(dims[1]) if box_ok else "",
            "box_z_A": float(dims[2]) if box_ok else "",
            "box_alpha_deg": float(dims[3]) if box_ok else "",
            "box_beta_deg": float(dims[4]) if box_ok else "",
            "box_gamma_deg": float(dims[5]) if box_ok else "",
        })
    return summary_rows, contact_rows


def write_csv(path: str | Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    """Write rows to CSV, creating the parent directory if needed."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
