# Protein Neighbor Analysis

**Version 0.1.0** — reusable neighborhood, water-contact, and residue-occupancy analysis for molecular-dynamics trajectories.

This repository refactors the useful methods from the historical `neighbor_test.zip` research workspace into a configurable Python package. It supports both direct analysis of MD trajectories and aggregation of the old frame-keyed JSON files, without depending on hardcoded protein names, residue numbers, phase names, or file paths.

## What it measures

For every analyzed frame, the trajectory workflow exports four distinct quantities:

- **Atom-pair contacts**: every central-selection/search-selection atom pair within the cutoff.
- **Unique neighboring atoms**: distinct atom indices participating in at least one contact.
- **Unique neighboring residues**: distinct topology residues participating in at least one contact. The internal residue index is included so residues with duplicate printed labels are not merged.
- **Frame/time/box metadata**: original frame index, trajectory time as reported by MDAnalysis, an explicitly labelled time unit, and whether periodic boundary conditions (PBC) were used.

Distances and cutoffs are in **Å** (the standard MDAnalysis coordinate unit). PBC defaults to `auto`: use the minimum-image convention only when a valid unit cell is present. Set `--pbc on` to require box dimensions or `--pbc off` for direct Cartesian distance.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate                 # Windows: .venv\Scripts\activate
python -m pip install -e ".[all,test]"
pytest -q
pna --help
```

Installations can be smaller when needed: `pip install -e ".[trajectory]"` for MDAnalysis trajectory processing, or `pip install -e ".[plot]"` for plots of archived JSON results. Core JSON summarization requires only NumPy and the Python standard library.

## Quick start: trajectory analysis

The repository contains a tiny synthetic PDB/DCD example (no research trajectory is bundled). Run from the repository root:

```bash
pna trajectory \
  --topology examples/toy_system.pdb \
  --trajectory examples/toy_trajectory.dcd \
  --center-selection "segid PROT and resid 1 and name SG" \
  --search-selection "not (segid PROT and resid 1)" \
  --cutoff 3.5 \
  --pbc auto \
  --output-dir neighbor_output \
  --prefix toy
```

Outputs:

- `neighbor_output/toy_per_frame.csv` — one row per analyzed frame.
- `neighbor_output/toy_contacts.csv` — one row per contacting atom pair, including center and neighbor atom identities and pair distance.
- `neighbor_output/toy_metadata.json` — run settings, selections, cutoff, frame range, and input filenames for reproducibility.

Add `--json` to export the detailed contact rows in JSON form. Change the selections, cutoff, frame range (`--start`, `--stop`, `--step`), and PBC mode to suit the system. Atom selections use the MDAnalysis selection language; see the [selection documentation](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).

The `--time-unit` setting labels the `ts.time` values supplied by the trajectory reader; it does **not** rescale them. Verify the trajectory time step and any external offset before interpreting time in physical units.

## Summarize the archived neighbor JSON results

The legacy files in `neighbor_test.zip` use the shape `{"0000": [{"segid": "...", "resid": 12, "resname": "CYS", "name": "O"}, ...], ...}`. To summarize one or several of these files:

```bash
pna summarize-json \
  --inputs path/to/phase_1.json path/to/phase_2.json path/to/phase_3.json \
  --output-dir json_summary \
  --plots
```

To count only selected neighbor atom names (following the spirit of the historical “electronegative atoms” plots), add a filter such as `--atom-names O,N,SG,ND1,ND2,OD1,OD2,OE1,OE2,OG,OG1,NE2,SD`. This is an **atom-name filter**, not a validated chemical/electronegativity classification; atom naming is force-field/topology dependent.

Each input file is kept as a separate series. The outputs are `frame_summary.csv`, `atom_occupancy.csv`, and `residue_occupancy.csv`, plus (if `--plots` is supplied) per-source neighbor-count traces, neighbor-count histograms, and residue-occupancy plots. The parser retains both source filename and source-frame key rather than inventing continuity between phases. Repeated basenames are automatically disambiguated with an input number, but distinct filenames are clearer.

## Scientific and data caveats

1. A contact cutoff is an analysis choice, not a universal physical boundary. Inspect sensitivity to plausible cutoffs and report it.
2. A central selection can contain multiple atoms. The code reports atom-pair count separately from the number of unique neighboring atoms/residues.
3. `--search-selection` is user-controlled. For water-only neighbors, select water atoms explicitly using residue/atom names that match the topology; do not assume all non-protein atoms are water (that can include ions, cofactors, and ligands).
4. PBC correctness depends on the trajectory containing valid box dimensions for each frame. The code records whether a box was actually used.
5. Legacy JSON records do not carry a unique topology `resindex` or necessarily a chain/insertion code, so residue identities in those files are best-effort based on segment, chain, printed residue id, and residue name.
6. The archived JSON/PNG/PDF products are prior analysis outputs, not a validated benchmark set. They are not included en masse. No original topology/trajectory inputs were included in the source archive, so this release is tested against synthetic trajectories, not a reproduction of the old published figures.
7. “Electronegative neighbors” should not be interpreted as a formal chemistry classification solely from an atom-name list. Prefer selection by topology element metadata when that metadata is present and trustworthy.

## Layout

```text
src/protein_neighbor_analysis/
  trajectory.py      # PBC-aware neighbor/contact algorithm and CSV writer
  legacy_json.py     # legacy JSON import, per-frame counts, residue occupancy
  plotting.py        # optional plots
  cli.py             # pna command line
examples/             # tiny synthetic trajectory and legacy JSON example
tests/                # unit and end-to-end checks
docs/                 # source audit, method definitions, validation report
```

## Validation

Run `pytest -q`. Tests verify PBC on/off behavior, invalid selections and cutoffs, empty selections, duplicate residue labels, correct source separation when summarizing multiple phase JSONs, atom-name filtering, and CLI output generation. The toy DCD/PDB provides an end-to-end example without exposing research data.

## Citation and license

If this workflow supports scientific work, cite MDAnalysis and your system-specific simulation/source paper. Recommended software references and reporting fields are listed in [`docs/references.md`](docs/references.md). The source archive contains research scripts and generated results authored for a particular project; please review the license/provenance note before public redistribution. No open-source license is declared in this release yet; see [`DISTRIBUTION-STATUS.md`](DISTRIBUTION-STATUS.md).
