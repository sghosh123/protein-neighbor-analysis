"""Command-line interface for protein-neighbor-analysis."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .legacy_json import summarize_legacy_json_files
from .trajectory import CONTACT_FIELDS, SUMMARY_FIELDS, analyze_universe, write_csv


def _write_contacts_json(path: Path, rows: list[dict], metadata: dict) -> None:
    payload = {
        "schema_version": 1,
        "distance_unit": "angstrom",
        "metadata": metadata,
        "contacts": rows,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _trajectory_command(args: argparse.Namespace) -> int:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError("MDAnalysis is needed for trajectory mode. Install with pip install '.[trajectory]'") from exc
    top = Path(args.topology)
    traj = Path(args.trajectory)
    if not top.is_file() or not traj.is_file():
        raise FileNotFoundError(f"Topology or trajectory missing: topology={top}, trajectory={traj}")
    universe = mda.Universe(str(top), str(traj))
    summary, contacts = analyze_universe(
        universe, args.center_selection, args.search_selection, args.cutoff,
        start=args.start, stop=args.stop, step=args.step, pbc=args.pbc,
        exclude_center_atoms=not args.include_center_atoms, time_unit=args.time_unit,
    )
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix
    summary_path = out_dir / f"{prefix}_per_frame.csv"
    contacts_path = out_dir / f"{prefix}_contacts.csv"
    metadata_path = out_dir / f"{prefix}_metadata.json"
    metadata = {
        "schema_version": 1,
        "topology_filename": top.name,
        "trajectory_filename": traj.name,
        "center_selection": args.center_selection,
        "search_selection": args.search_selection,
        "cutoff_A": float(args.cutoff),
        "pbc_requested": args.pbc,
        "exclude_center_atoms_from_search": not args.include_center_atoms,
        "start_frame_inclusive": args.start,
        "stop_frame_exclusive": args.stop,
        "frame_step": args.step,
        "time_unit_label": args.time_unit,
        "time_values_note": "Copied from MDAnalysis Timestep.time; not rescaled by this tool.",
        "frames_available": int(universe.trajectory.n_frames),
        "frames_analyzed": len(summary),
    }
    write_csv(summary_path, summary, SUMMARY_FIELDS)
    write_csv(contacts_path, contacts, CONTACT_FIELDS)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if args.json:
        _write_contacts_json(out_dir / f"{prefix}_contacts.json", contacts, metadata)
    print(f"Frames analyzed: {len(summary)}")
    print(f"Per-frame metrics: {summary_path}")
    print(f"Atom-pair contacts: {contacts_path} ({len(contacts)} rows)")
    print(f"Run metadata: {metadata_path}")
    if args.json:
        print(f"Contact JSON: {out_dir / f'{prefix}_contacts.json'}")
    return 0


def _json_command(args: argparse.Namespace) -> int:
    atom_names = None
    if args.atom_names:
        atom_names = {x.strip().upper() for x in args.atom_names.split(",") if x.strip()}
        if not atom_names:
            raise ValueError("--atom-names was provided but contains no usable names")
    frame_rows, atom_rows, occupancy_rows = summarize_legacy_json_files(args.inputs, atom_names=atom_names)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame_path = out / "frame_summary.csv"
    atom_path = out / "atom_occupancy.csv"
    occupancy_path = out / "residue_occupancy.csv"
    with frame_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(frame_rows[0].keys()) if frame_rows else [
            "source_file", "source_frame", "frame_ordinal_within_file", "n_atom_records",
            "n_unique_neighbor_atoms", "n_unique_neighbor_residues", "filter_atom_names"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(frame_rows)
    with atom_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["source_file", "segid", "chain_id", "resid", "resname", "atom_name", "frames_present", "total_frames", "occupancy_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(atom_rows)
    with occupancy_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["source_file", "segid", "chain_id", "resid", "resname", "frames_present", "total_frames", "occupancy_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(occupancy_rows)
    print(f"Frame summary: {frame_path} ({len(frame_rows)} frames)")
    print(f"Atom occupancy: {atom_path} ({len(atom_rows)} source/atom rows)")
    print(f"Residue occupancy: {occupancy_path} ({len(occupancy_rows)} source/residue rows)")
    if args.plots:
        from .plotting import plot_json_summaries
        for path in plot_json_summaries(frame_rows, occupancy_rows, out / "plots"):
            print(f"Plot: {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pna", description="Protein/atom neighbor analysis for MD trajectories and archived JSON results")
    sub = parser.add_subparsers(dest="command", required=True)
    traj = sub.add_parser("trajectory", help="Calculate neighbors directly from an MDAnalysis topology/trajectory")
    traj.add_argument("--topology", required=True, help="Topology file understood by MDAnalysis (for example PSF or PDB)")
    traj.add_argument("--trajectory", required=True, help="Trajectory file understood by MDAnalysis (for example DCD or XTC)")
    traj.add_argument("--center-selection", required=True, help="MDAnalysis selection defining the central atom(s)")
    traj.add_argument("--search-selection", default="all", help="MDAnalysis selection defining candidates to search (default: all)")
    traj.add_argument("--cutoff", required=True, type=float, help="Distance cutoff in Angstrom")
    traj.add_argument("--start", type=int, default=0, help="First frame index (inclusive; default: 0)")
    traj.add_argument("--stop", type=int, default=None, help="Frame stop index (exclusive; default: end)")
    traj.add_argument("--step", type=int, default=1, help="Analyze every Nth frame (default: 1)")
    traj.add_argument("--pbc", choices=["auto", "on", "off"], default="auto", help="Periodic boundary handling (default: auto)")
    traj.add_argument("--time-unit", default="ps", help="Label for the trajectory time values (default: ps; does not rescale times)")
    traj.add_argument("--include-center-atoms", action="store_true", help="Do not exclude central atom indices from the search group")
    traj.add_argument("--output-dir", default="neighbor_output", help="Output directory")
    traj.add_argument("--prefix", default="neighbors", help="Output filename prefix")
    traj.add_argument("--json", action="store_true", help="Also export flat atom-pair contact records as JSON")
    traj.set_defaults(func=_trajectory_command)

    legacy = sub.add_parser("summarize-json", help="Aggregate legacy frame -> atom-record neighbor JSON files")
    legacy.add_argument("--inputs", nargs="+", required=True, help="One or more JSON files; each remains an independent series")
    legacy.add_argument("--output-dir", default="neighbor_json_summary", help="Output directory")
    legacy.add_argument("--atom-names", default=None, help="Optional comma-separated atom-name filter (e.g. O,N,SG,OD1)")
    legacy.add_argument("--plots", action="store_true", help="Write neighbor-count traces and top residue-occupancy plots")
    legacy.set_defaults(func=_json_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, OSError, RuntimeError, ImportError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
