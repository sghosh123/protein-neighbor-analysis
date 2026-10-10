# Release validation — v0.1.0

## Test environment

- Python 3.13.5
- MDAnalysis 2.10.0
- NumPy 2.3.5
- Matplotlib 3.10.8

## Checks performed

- `pytest -q`: **16 passed**.
- `python -m compileall -q src tests`: passed.
- Editable installation with the declared extras: passed.
- Non-editable Python wheel build: passed.
- Wheel installed to an isolated target directory and `python -m protein_neighbor_analysis --help` executed successfully.
- Wheel contents checked to ensure the distribution note was not packaged as if it were a software license.
- Installed CLI examples executed against the synthetic PDB/DCD and legacy JSON fixtures. Trajectory run: 3 frames, 6 atom-pair contact records. Legacy JSON run: 3 frames, five atom-occupancy rows, four residue-occupancy rows, and three PNG plots.

The release ZIP is checked separately for CRC integrity, unwanted cache/build files, and presence of both synthetic fixtures. A sibling `.sha256` file is provided for the exact downloadable archive.
