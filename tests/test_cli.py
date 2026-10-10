import csv
import json
from pathlib import Path

import pytest
from protein_neighbor_analysis.cli import main


def test_summarize_json_cli_outputs_tables_and_plots(tmp_path):
    source = tmp_path / "neighbors.json"
    source.write_text(json.dumps({
        "0": [{"segid": "PROT", "resid": 2, "resname": "ALA", "name": "O"}],
        "1": [],
    }), encoding="utf-8")
    out = tmp_path / "out"
    assert main(["summarize-json", "--inputs", str(source), "--output-dir", str(out), "--plots"]) == 0
    frame_path = out / "frame_summary.csv"
    occ_path = out / "residue_occupancy.csv"
    assert frame_path.exists() and occ_path.exists()
    assert (out / "atom_occupancy.csv").exists()
    with frame_path.open(newline="", encoding="utf-8") as h:
        rows = list(csv.DictReader(h))
    assert len(rows) == 2
    assert rows[0]["n_unique_neighbor_atoms"] == "1"
    assert len(list((out / "plots").glob("*.png"))) >= 3


@pytest.mark.filterwarnings("ignore:DCDReader currently makes independent timesteps:DeprecationWarning")
def test_trajectory_cli_emits_summary_and_contacts(tmp_path):
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    top = tmp_path / "toy.pdb"
    traj = tmp_path / "toy.dcd"
    # Create the files with an in-memory topology; positions are separated by 1 Å.
    u = mda.Universe.empty(2, n_residues=2, atom_resindex=np.arange(2), trajectory=True)
    u.add_TopologyAttr("names", ["SG", "O"])
    u.add_TopologyAttr("types", ["S", "O"])
    u.add_TopologyAttr("elements", ["S", "O"])
    u.add_TopologyAttr("resnames", ["CYS", "TIP3"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("segids", ["SYS"])
    u.add_TopologyAttr("ids", [1, 2])
    u.add_TopologyAttr("altLocs", [" ", " "])
    u.add_TopologyAttr("icodes", [" ", " "])
    u.add_TopologyAttr("chainIDs", ["A", "A"])
    u.add_TopologyAttr("occupancies", [1.0, 1.0])
    u.add_TopologyAttr("tempfactors", [0.0, 0.0])
    u.add_TopologyAttr("record_types", ["ATOM", "HETATM"])
    u.add_TopologyAttr("formalcharges", [0, 0])
    coords = np.array([[[5, 5, 5], [5.7, 5, 5]], [[5, 5, 5], [7, 5, 5]]], dtype=np.float32)
    u.load_new(coords, format=MemoryReader, dimensions=np.tile([10, 10, 10, 90, 90, 90], (2, 1)))
    with mda.coordinates.PDB.PDBWriter(str(top), multiframe=False, reindex=True) as writer:
        writer.write(u.atoms)
    with mda.Writer(str(traj), n_atoms=2) as writer:
        for ts in u.trajectory:
            writer.write(u.atoms)
    out = tmp_path / "analysis"
    assert main(["trajectory", "--topology", str(top), "--trajectory", str(traj),
                 "--center-selection", "name SG", "--search-selection", "name O",
                 "--cutoff", "1.0", "--pbc", "off", "--output-dir", str(out)]) == 0
    assert (out / "neighbors_per_frame.csv").exists()
    assert (out / "neighbors_contacts.csv").exists()
    metadata = json.loads((out / "neighbors_metadata.json").read_text(encoding="utf-8"))
    assert metadata["cutoff_A"] == 1.0
    assert metadata["center_selection"] == "name SG"
    assert metadata["frames_analyzed"] == 2
