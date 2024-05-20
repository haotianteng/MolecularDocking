#run rdkit on pdbbind ligand as a test
import os
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from mdbm.utils.mol import read_molecule

def run_optimize(f):
    mol = read_molecule(f,sanitize=False, remove_hs=False, calc_charges=False)
    if mol is None:
        raise ValueError(f"Failed to read {f}")
    mol = Chem.AddHs(mol)
    # AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    out_f = f.replace(".mol2", "_rdkit.sdf")
    w = Chem.SDWriter(out_f)
    w.write(mol)

if __name__ == "__main__":
    df = 'datasets/PDBBind_processed_test/'
    pdb_ids = os.listdir(df)
    for pdb_id in tqdm(pdb_ids):
        sub_f = os.path.join(df, pdb_id)
        ligand_f = os.path.join(sub_f, f'{pdb_id}_ligand.mol2')
        run_optimize(ligand_f)