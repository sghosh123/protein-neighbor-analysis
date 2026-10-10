# Integration report — protein-neighbor-analysis v0.1.0

## Source basis

Primary source: `neighbor_test.zip` (SHA-256 `86a9dfea18135c5d5eb2458aaafdfd89688d46bc45986e507996f34096a5d95e`). The archive was 450 entries (~110 MB uncompressed) and contained nine Python scripts, 101 frame-keyed neighbor JSON files, 339 PNGs, and one PDF. The archive passed the ZIP integrity check.

## Integrated workflows

- PBC-aware per-frame neighbor analysis using explicit MDAnalysis center/search selections.
- Atom-pair, unique-atom, and unique-residue contact metrics in separate output fields.
- CSV contact detail with central/neighbor identities, distance in Å, topology residue index, time labels, and frame metadata; a run-metadata JSON captures selections, cutoffs, frame sampling, and input filenames.
- Legacy JSON import and validation, per-frame count metrics, atom occupancy, residue occupancy, optional atom-name filtering, and per-source trace/histogram/residue-occupancy plots.
- Synthetic PDB/DCD example and a small legacy JSON fixture; no original research trajectories or large historical result exports bundled.

## Corrections incorporated

- Fixed the hardcoded center/search selection behavior observed in `find_neighbor_waters.py` by using actual user-provided selections.
- Prevented phase-data reuse in legacy aggregation by reading and preserving every input dataset independently.
- Distinguished atom-pair count, unique atom count, and unique residue count; topology residue indices prevent duplicate human-readable labels from merging in trajectory mode.
- Added explicit `auto`, `on`, and `off` PBC modes and records whether a valid box was applied for each analyzed frame.
- Replaced fragile histogram-bin calculations based on `np.ptp` with integer-centered bins that work for constant-valued counts.
- Made atom-name filtering case-insensitive and disambiguated duplicate input basenames.
- Documented limits in old JSON identity fields, cutoffs, time units, water selections, and the “electronegative atom” label.
- Did not copy the legacy script with mixed indentation that fails Python parsing.

## Validation performed

- `pytest -q`: 16 passed.
- `python -m compileall -q src tests`: passed.
- Editable package installation: passed.
- Standard wheel build: passed.
- Installed CLI help: passed.
- End-to-end toy trajectory command: 3 frames processed, per-frame CSV and contact CSV/JSON written (6 atom-pair contact rows for the documented 3.5 Å example).
- End-to-end legacy JSON command: 3 frames summarized; frame, atom-occupancy, residue-occupancy CSVs and 3 plot files written.
- Release archive is distributed with an external SHA-256 sidecar; see `docs/release-validation.md` for the tested version/environment details.

Tests use synthetic data; no historical research result has been claimed reproduced. Before scientific use, validate selections, distance units, PBC dimensions, saved-frame spacing, and sample frame distances with the intended topology/trajectory pair.
