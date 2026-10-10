# Workflow catalog

| Workflow | Entry point | Output | Dependencies |
|---|---|---|---|
| Trajectory neighbor analysis | `pna trajectory` | Per-frame CSV + atom-pair CSV + run metadata JSON (+ optional contact JSON) | MDAnalysis |
| Archive JSON summary | `pna summarize-json` | Per-frame counts + residue occupancy CSV | Standard library; NumPy installed with package |
| Archive summary plots | `pna summarize-json --plots` | Per-source traces and occupancy plots | matplotlib |

## Important user parameters

`center_selection` defines the geometric center atom set (a single atom or group); `search_selection` defines candidates. `cutoff` is Å. `pbc=auto` uses a valid box on each frame when available; `on` requires valid box data; `off` intentionally disables PBC. `start`/`stop`/`step` follow Python slicing semantics: start inclusive, stop exclusive.

The result summarizes the selected atom set's geometry only. It does not infer binding, hydrogen bonds, energetic favorability, or chemical coordination by itself.
