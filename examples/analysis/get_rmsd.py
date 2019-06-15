import tempfile
import pandas as pd
import sys
import pickle
import os
import glob
import multiprocessing
import cpeptools
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import mdtraj as md



def get_rmsd(smiles, pdb_path, ref_pdb_path):
    smiles_mol = Chem.MolFromSmiles(smiles)
    ref_mol = Chem.MolFromPDBFile(ref_pdb_path)
    mol = Chem.MolFromPDBFile(pdb_path)

    ref_mol = AllChem.AssignBondOrdersFromTemplate(smiles_mol,ref_mol)
    mol = AllChem.AssignBondOrdersFromTemplate(smiles_mol, mol)
    order = list(mol.GetSubstructMatches(ref_mol)[0])
    mol = Chem.RenumberAtoms(mol, order)

    indices = cpeptools.get_largest_ring(ref_mol)
    assert len(set(indices) - set(cpeptools.get_largest_ring(mol))) == 0, "ring atom indices do not agree"

    tmp_dir = tempfile.mkdtemp()
    ref_pdb_filename = tempfile.mktemp(suffix=".pdb", dir = tmp_dir)
    pdb_filename = tempfile.mktemp(suffix=".pdb", dir = tmp_dir)
    Chem.MolToPDBFile(ref_mol, ref_pdb_filename)
    Chem.MolToPDBFile(mol, pdb_filename)

    ref  = md.load(ref_pdb_filename)
    compare = md.load(pdb_filename)


    rmsd = md.rmsd(compare, ref, 0)
    ring_rmsd = md.rmsd(compare, ref, 0, atom_indices = indices)
    compare = compare.superpose(ref, 0, atom_indices = indices)
    return rmsd, compare[np.argmin(rmsd)] , ring_rmsd, compare[np.argmin(ring_rmsd)] , ref


import glob
def calculate(name):
    if len(glob.glob("../conformer/{}/{}/*.pdb".format(variant, name))) == 0:
        print("{}/{} has no conformers".format(variant, name))
        return
    if not os.path.exists(variant + "/" + name):
        os.makedirs(variant + "/" + name)
    confs = md.load(glob.glob("../conformer/{}/{}/*pdb".format(variant,name)))
    print(" {} has {} conformers".format(name, len(confs)))

    try:
        tmp_dir = tempfile.mkdtemp()
        tmp_pdb_filename = tempfile.mktemp(suffix=".pdb", dir = tmp_dir)
        confs.save(tmp_pdb_filename)

        print("starting {}".format(name))
        rmsd, rmsd_pdb, ring_rmsd, ring_rmsd_pdb, ref =  get_rmsd(mol_dict[name], tmp_pdb_filename, '../reference/{}.pdb'.format(name))

        ring_rmsd_pdb.save("{}/{}/ring_rmsd_min.pdb".format(variant, name))
        np.save("{}/{}/ring_rmsd.npy".format(variant, name), ring_rmsd)
        rmsd_pdb.save("{}/{}/rmsd_min.pdb".format(variant, name))
        np.save("{}/{}/rmsd.npy".format(variant, name), rmsd)
        ref.save("{}/{}/aligned_ref.pdb".format(variant, name))
    except Exception as e:
        print("Problem with {} : {}".format(name, e))


# the folder name in conformer, whether it is mETKDG, mETKDG with
variant = sys.argv[1]
names = list(os.walk("../conformer/{}/".format(variant)))[0][1]
df = pd.read_csv("../reference/smiles.tsv", sep = "\t", comment = "%")

mol_dict = {name : val for name, val in zip(df.Name, df.Smiles)}

for name in names:
    try:
        calculate(name)
    except:
        print("Problem with {}".format(name))
