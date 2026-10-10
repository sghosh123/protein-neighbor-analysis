# Source audit: `neighbor_test.zip`

## Source

- Archive: `neighbor_test.zip`
- SHA-256: `86a9dfea18135c5d5eb2458aaafdfd89688d46bc45986e507996f34096a5d95e`
- Archive validation: `zipfile.testzip()` reported no corrupt entry.
- Inventory: 450 entries, including 9 Python scripts, 101 JSON datasets, 339 PNGs, and one PDF (~110 MB uncompressed).

The release intentionally does not copy the full result bundle. It keeps the generalizable method as code and uses synthetic fixtures so that the repository remains small and independent of unpublished system-specific data.

## Bugs and assumptions identified

1. **Ignored function arguments:** `yeast_MBD_apo_neighbor/find_neighbor_waters.py` builds the requested central/search selection strings but actually calls `select_atoms("not protein")` and `select_atoms("protein and resid 16 and name HG1")`. Thus `protein_segname` and `protein_resid` do not control the calculation. The new CLI uses the actual center and search selection arguments.
2. **Phase data reused accidentally:** `binding_loop_cys_neighbor_list/human/electronegative_neighbor_count_by_residue.py` loads four phase files but builds `value2`, `value3`, and `value4` from `dictionarya`, not dictionaries b/c/d. This can give invalid phase-combined contact trends. The new JSON summarizer reads each input file independently and retains source identity.
3. **Syntax error:** `yeast_MBD_apo_neighbor/water_neighbor_count_by_residue.py` uses inconsistent tabs/spaces and fails Python parsing under standard Python 3.
4. **Ambiguous counts:** legacy plots variously use raw atom entries, de-duplicated atom names, residue numbers, or named atom/residue strings. Those quantities were all described as “contacts” in some plots. The new output names atom-pair contacts, unique neighbor atoms, and unique neighbor residues separately.
5. **Hardcoded scientific settings:** old scripts embed filenames, protein/residue ids, phase count, cutoffs, titles, plot limits, and frame extents in code. Histogram bin counts often use `np.ptp(counts)`, which becomes zero for constant-valued data; plots sometimes use an axes object from another figure.
6. **Potential residue identity collisions:** legacy JSON lacks a topology `resindex`, insertion code, or guaranteed chain field. A JSON aggregation can only distinguish residues as far as the exported metadata permits. New trajectory analysis uses `resindex` as an unambiguous within-topology identifier.
7. **Selection is not a chemistry definition:** `not protein` can include ions, cofactors, solvent, and other non-protein atoms. Water-only analysis must select the water residue/atom names used in the supplied topology.
8. **Atom-name “electronegative” filters:** the original explicit atom-name list is force-field/topology dependent and is not equivalent to a general element classification. The CLI describes it as an optional atom-name filter.
9. **PBC and units:** historical methods pass per-frame dimensions to MDAnalysis neighbor search. The new method validates unit-cell dimensions on each frame, supports `auto/on/off`, and records whether a box was applied. All cutoffs/distances are Å; time units are only labels for trajectory-reader time values and are never silently scaled.

## What was retained

- The per-frame central atom/residue neighborhood idea.
- Optional filtering to selected neighbor atom names.
- Neighbor-count traces, histograms, and residue occupancy/population summaries.
- Legacy JSON import so previous results can be examined without recomputing trajectories.

## What was not carried over as-is

The 101 JSON data products, 339 plot images, and large PDF are not copied into the reusable package. Their research provenance, source topology/trajectory pairing, and publication/distribution permissions should be checked before any result bundle is released. No matching raw trajectories or topologies were part of this archive, so the original results cannot be reproduced from it alone.
