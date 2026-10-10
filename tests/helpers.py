import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader


def make_universe(coords=None, dimensions=None, *, duplicate_residue_labels=False):
    """Build a small in-memory topology and trajectory for deterministic tests."""
    if coords is None:
        coords = np.array([[[9.8, 5.0, 5.0], [0.4, 5.0, 5.0], [0.6, 5.0, 5.0]]], dtype=np.float32)
    coords = np.asarray(coords, dtype=np.float32)
    n_atoms = coords.shape[1]
    if duplicate_residue_labels:
        atom_resindex = np.arange(n_atoms)
        n_residues = n_atoms
        resids = [1, 2, 2][:n_atoms]
        resnames = ["CYS", "ALA", "ALA"][:n_atoms]
        names = ["SG", "O", "O"][:n_atoms]
        elements = ["S", "O", "O"][:n_atoms]
    else:
        atom_resindex = np.arange(n_atoms)
        n_residues = n_atoms
        resids = list(range(1, n_atoms + 1))
        resnames = ["CYS"] + ["TIP3"] * (n_atoms - 1)
        names = ["SG"] + ["O"] * (n_atoms - 1)
        elements = ["S"] + ["O"] * (n_atoms - 1)
    u = mda.Universe.empty(n_atoms, n_residues=n_residues, atom_resindex=atom_resindex, trajectory=True)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("types", names)
    u.add_TopologyAttr("elements", elements)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", resids)
    u.add_TopologyAttr("segids", ["SYS"])
    if dimensions is None:
        dimensions = np.tile(np.array([10, 10, 10, 90, 90, 90], dtype=np.float32), (coords.shape[0], 1))
    u.load_new(coords, format=MemoryReader, dimensions=np.asarray(dimensions, dtype=np.float32))
    return u
