"""Summaries and exports for legacy per-frame neighbor JSON data."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


def _atom_key(atom: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        str(atom.get("segid", "")), str(atom.get("chain_id", atom.get("chainID", ""))),
        str(atom.get("resid", "")), str(atom.get("resname", "")), str(atom.get("name", "")),
    )


def _residue_key(atom: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(atom.get("segid", "")), str(atom.get("chain_id", atom.get("chainID", ""))),
        str(atom.get("resid", "")), str(atom.get("resname", "")),
    )


def read_legacy_json(path: str | Path) -> dict[str, list[dict[str, Any]]]:
    """Read legacy frame -> list-of-atom-records JSON and validate its shape."""
    source = Path(path)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read JSON {source}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object keyed by frame in {source}")
    cleaned: dict[str, list[dict[str, Any]]] = {}
    for frame, neighbors in payload.items():
        if not isinstance(neighbors, list) or any(not isinstance(x, dict) for x in neighbors):
            raise ValueError(f"Frame {frame!r} in {source} must contain a list of atom objects")
        cleaned[str(frame)] = neighbors
    return cleaned


def summarize_legacy_json_files(
    paths: Iterable[str | Path], *, atom_names: set[str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Summarize each input JSON independently and calculate residue occupancy.

    Input files remain separate series in the output. This prevents the
    historical mistake of accidentally reusing phase-1 data for later phases.
    """
    files = [Path(p) for p in paths]
    if not files:
        raise ValueError("Provide at least one input JSON file")
    if atom_names is not None:
        atom_names = {str(name).strip().upper() for name in atom_names if str(name).strip()}
        if not atom_names:
            raise ValueError("atom_names must contain at least one non-empty atom name")
    frame_rows: list[dict[str, Any]] = []
    atom_occupancy: dict[tuple[str, tuple[str, str, str, str, str]], set[str]] = defaultdict(set)
    residue_occupancy: dict[tuple[str, tuple[str, str, str, str]], set[str]] = defaultdict(set)
    totals: dict[str, int] = defaultdict(int)
    used_names: set[str] = set()
    for file_number, file_path in enumerate(files, start=1):
        data = read_legacy_json(file_path)
        source_name = file_path.name
        if source_name in used_names:
            source_name = f"{file_path.stem}__input{file_number}{file_path.suffix}"
            while source_name in used_names:
                source_name = f"{file_path.stem}__input{file_number}_{len(used_names)}{file_path.suffix}"
        used_names.add(source_name)
        totals[source_name] = len(data)
        for ordinal, (source_frame, atoms) in enumerate(data.items()):
            filtered = atoms if atom_names is None else [a for a in atoms if str(a.get("name", "")).upper() in atom_names]
            unique_atoms = {_atom_key(a) for a in filtered}
            unique_residues = {_residue_key(a) for a in filtered}
            frame_rows.append({
                "source_file": source_name,
                "source_frame": str(source_frame),
                "frame_ordinal_within_file": ordinal,
                "n_atom_records": len(filtered),
                "n_unique_neighbor_atoms": len(unique_atoms),
                "n_unique_neighbor_residues": len(unique_residues),
                "filter_atom_names": "" if atom_names is None else ";".join(sorted(atom_names)),
            })
            for atom_key in unique_atoms:
                atom_occupancy[(source_name, atom_key)].add(str(source_frame))
            for residue in unique_residues:
                residue_occupancy[(source_name, residue)].add(str(source_frame))

    atom_rows: list[dict[str, Any]] = []
    for (source_name, atom_key), frames in sorted(atom_occupancy.items()):
        total = totals[source_name]
        segid, chain_id, resid, resname, atom_name = atom_key
        atom_rows.append({
            "source_file": source_name, "segid": segid, "chain_id": chain_id,
            "resid": resid, "resname": resname, "atom_name": atom_name,
            "frames_present": len(frames), "total_frames": total,
            "occupancy_fraction": len(frames) / total if total else 0.0,
        })

    residue_rows: list[dict[str, Any]] = []
    for (source_name, residue), frames in sorted(residue_occupancy.items()):
        total = totals[source_name]
        segid, chain_id, resid, resname = residue
        residue_rows.append({
            "source_file": source_name, "segid": segid, "chain_id": chain_id,
            "resid": resid, "resname": resname, "frames_present": len(frames),
            "total_frames": total, "occupancy_fraction": len(frames) / total if total else 0.0,
        })
    return frame_rows, atom_rows, residue_rows
