# Methodology and output definitions

## Geometric contact search

For each analyzed frame, the central and candidate atom selections are resolved once from the topology. The implementation then computes distances between the selected coordinate arrays using `MDAnalysis.lib.distances.capped_distance`. With a valid unit cell and PBC enabled, MDAnalysis applies the minimum-image convention; otherwise `pbc=auto` records and uses direct distances. `pbc=on` fails if a frame has no valid box, rather than silently running without PBC.

A pair is considered a contact when its distance is less than or equal to the configured maximum cutoff, as implemented by MDAnalysis. The cutoff is in Å. The frame metrics are:

- `n_contact_pairs`: number of center/candidate atom pairs inside cutoff.
- `n_neighbor_atoms`: number of distinct candidate atom indices inside cutoff of at least one center atom.
- `n_neighbor_residues`: number of distinct topology residue indices represented among neighbor atoms.

These are intentionally different numbers. In particular, one central atom can contact several atoms in the same water or residue.

## Residue identities

For topology-backed trajectory analysis, the unique MDAnalysis `resindex` is part of the identity key. This prevents two distinct topology residues with duplicated `(segid, resid, resname)` labels from being merged. Human-readable labels are also exported.

Legacy neighbor JSON files were written without the unique `resindex`, and often omit chain identifiers and insertion codes. Occupancy derived from them necessarily falls back to the best available printed identifiers; the provenance limitation is documented in the report.

## Time metadata

`time` is copied from `ts.time`; `--time-unit` only sets the label. It does not change the underlying value or apply any time offset. Trajectory formats and the metadata used when loading them determine the reported values. Check `dt`, offsets, saved-frame spacing, and trajectory origin before interpreting a count as a duration in ns.

## Water and electronegative atom analysis

The package does not guess that `not protein` means water. Users should supply a search selection aligned with their topology (for example `resname TIP3 and name OH2`, if those are the actual names). The archived atom-name list for “electronegative” neighbors is exposed as an optional filter, not treated as a universal chemistry rule; element-based classification is preferable when topology element metadata is valid.
