import json

import pytest

from protein_neighbor_analysis.legacy_json import read_legacy_json, summarize_legacy_json_files


def _write(path, data):
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_multiple_phase_files_are_summarized_independently(tmp_path):
    p1 = _write(tmp_path / "phase1.json", {"0000": [{"segid": "A", "resid": 1, "resname": "CYS", "name": "SG"}], "0001": []})
    p2 = _write(tmp_path / "phase2.json", {"0100": [{"segid": "A", "resid": 2, "resname": "HOH", "name": "O"}, {"segid": "A", "resid": 2, "resname": "HOH", "name": "H1"}]})
    frames, atom_occ, occupancy = summarize_legacy_json_files([p1, p2])
    assert len(frames) == 3
    assert frames[0]["n_unique_neighbor_atoms"] == 1
    assert frames[1]["n_unique_neighbor_atoms"] == 0
    assert frames[2]["n_unique_neighbor_atoms"] == 2
    assert {row["source_file"] for row in frames} == {"phase1.json", "phase2.json"}
    phase2_water = [r for r in occupancy if r["source_file"] == "phase2.json"][0]
    assert phase2_water["frames_present"] == 1
    assert phase2_water["total_frames"] == 1


def test_atom_name_filter_is_case_insensitive(tmp_path):
    path = _write(tmp_path / "neighbors.json", {"0": [
        {"segid": "WAT", "resid": 1, "resname": "TIP3", "name": "OH2"},
        {"segid": "WAT", "resid": 1, "resname": "TIP3", "name": "H1"},
    ]})
    frames, atom_occ, occupancy = summarize_legacy_json_files([path], atom_names={"oh2"})
    assert frames[0]["n_atom_records"] == 1
    assert frames[0]["n_unique_neighbor_residues"] == 1
    assert occupancy[0]["occupancy_fraction"] == 1.0
    assert atom_occ[0]["atom_name"] == "OH2"


def test_invalid_json_shape_is_rejected(tmp_path):
    p = _write(tmp_path / "bad.json", {"0": {"resid": 1}})
    with pytest.raises(ValueError, match="list of atom objects"):
        read_legacy_json(p)


def test_empty_frame_values_are_supported(tmp_path):
    p = _write(tmp_path / "empty.json", {"0": [], "1": []})
    frames, atom_occ, occupancy = summarize_legacy_json_files([p])
    assert len(frames) == 2
    assert all(r["n_unique_neighbor_atoms"] == 0 for r in frames)
    assert atom_occ == []
    assert occupancy == []


def test_duplicate_basenames_are_disambiguated(tmp_path):
    first_dir = tmp_path / "one"
    second_dir = tmp_path / "two"
    first_dir.mkdir(); second_dir.mkdir()
    p1 = _write(first_dir / "neighbors.json", {"0": []})
    p2 = _write(second_dir / "neighbors.json", {"0": [{"segid": "A", "resid": 1, "resname": "ALA", "name": "O"}]})
    frames, atom_occ, residue_occ = summarize_legacy_json_files([p1, p2])
    assert len({row["source_file"] for row in frames}) == 2
    assert frames[0]["source_file"] == "neighbors.json"
    assert frames[1]["source_file"].startswith("neighbors__input2")
