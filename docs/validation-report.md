# Validation report — version 0.1.0

The package is tested against controlled synthetic systems, not against the absent historical trajectories.

Coverage includes: valid/invalid cutoff values; selection failures; PBC minimum-image distances; direct non-PBC distances; missing box behavior; frame slicing; separation of atom-pair, atom, and residue counts; distinct topology residues with colliding printed residue labels; legacy JSON validation; empty frames; independent phase files; filtered atom-name counting; residue occupancy; CLI file generation including a run-metadata JSON; and plots from the sample JSON.

Before applying to research data, run with the exact topology and trajectory pair, verify MDAnalysis selections, inspect box dimensions and distance units, and check a few frame-level distances independently. No existing research result has been claimed reproduced.
