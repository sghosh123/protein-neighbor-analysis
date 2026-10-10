import numpy as np
import pytest

from protein_neighbor_analysis.trajectory import analyze_universe, valid_dimensions
from tests.helpers import make_universe


def test_pbc_wraps_across_box_boundary():
    u = make_universe()
    summary, contacts = analyze_universe(u, "name SG", "name O", 1.0, pbc="auto")
    assert summary[0]["pbc_used"] is True
    assert summary[0]["n_contact_pairs"] == 2
    assert summary[0]["n_neighbor_atoms"] == 2
    assert summary[0]["n_neighbor_residues"] == 2
    assert all(row["distance_A"] < 1.0 for row in contacts)


def test_nonperiodic_distance_does_not_wrap():
    u = make_universe()
    summary, contacts = analyze_universe(u, "name SG", "name O", 1.0, pbc="off")
    assert summary[0]["pbc_used"] is False
    assert summary[0]["n_contact_pairs"] == 0
    assert contacts == []


def test_duplicate_printed_residue_labels_are_not_merged():
    coords = np.array([[[5, 5, 5], [5.5, 5, 5], [5.7, 5, 5]]], dtype=np.float32)
    u = make_universe(coords, duplicate_residue_labels=True)
    summary, contacts = analyze_universe(u, "name SG", "name O", 1.0)
    assert summary[0]["n_neighbor_residues"] == 2
    assert len(contacts) == 2


def test_cutoff_must_be_positive_finite():
    u = make_universe()
    for cutoff in [0, -1, float("nan"), float("inf")]:
        with pytest.raises(ValueError, match="cutoff"):
            analyze_universe(u, "name SG", "name O", cutoff)


def test_empty_or_invalid_selection_is_reported():
    u = make_universe()
    with pytest.raises(ValueError, match="Center selection matched no atoms"):
        analyze_universe(u, "name ZZ", "name O", 2.0)
    with pytest.raises(ValueError, match="Invalid MDAnalysis selection"):
        analyze_universe(u, "name SG and (", "name O", 2.0)


def test_linear_unit_cell_angle_is_invalid():
    assert not valid_dimensions([10, 10, 10, 90, 90, 180])


def test_pbc_on_requires_a_valid_box():
    dims = np.zeros((1, 6), dtype=np.float32)
    u = make_universe(dimensions=dims)
    assert not valid_dimensions(u.trajectory.ts.dimensions)
    summary, _ = analyze_universe(u, "name SG", "name O", 2.0, pbc="auto")
    assert summary[0]["pbc_used"] is False
    with pytest.raises(ValueError, match="PBC was required"):
        analyze_universe(u, "name SG", "name O", 2.0, pbc="on")


def test_frame_slicing_preserves_original_frame_indices():
    coords = np.array([
        [[5, 5, 5], [5.5, 5, 5]],
        [[5, 5, 5], [7.0, 5, 5]],
        [[5, 5, 5], [5.7, 5, 5]],
    ], dtype=np.float32)
    u = make_universe(coords)
    summary, _ = analyze_universe(u, "name SG", "name O", 1.0, start=0, stop=3, step=2)
    assert [r["frame"] for r in summary] == [0, 2]
    assert [r["n_contact_pairs"] for r in summary] == [1, 1]


def test_multiple_center_atoms_count_atom_pairs_separately_from_unique_neighbors():
    coords = np.array([[[5.0, 5.0, 5.0], [5.3, 5.0, 5.0], [5.1, 5.0, 5.0]]], dtype=np.float32)
    u = make_universe(coords)
    u.atoms.names = ["SG1", "SG2", "O"]
    summary, contacts = analyze_universe(u, "name SG1 SG2", "name O", 0.5, pbc="off")
    assert summary[0]["n_contact_pairs"] == 2
    assert summary[0]["n_neighbor_atoms"] == 1
    assert summary[0]["n_neighbor_residues"] == 1
    assert len(contacts) == 2
