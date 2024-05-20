import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

def read_molecule(molecule_file, sanitize=False, calc_charges=False, remove_hs=False):
    if molecule_file.endswith('.mol2'):
        with open(molecule_file,'r') as file:
            mol = Chem.MolFromMol2Block(file.read(), sanitize=False, removeHs=False)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        raise ValueError('Expect the format of the molecule_file to be '
                         'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))

    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                print('Unable to compute charges for the molecule.')

        if remove_hs:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)

    except Exception as e:
        # Print stacktrace
        import traceback
        msg = traceback.format_exc()
        print(f"Failed to process molecule: {molecule_file}\n{msg}")
        return None
    
    return mol

def get_positions(mol):
    """
    Get the positions of the atoms in the molecule.
    """
    positions = []
    atoms = []
    conformer = mol.GetConformer()
    for i,atom in enumerate(mol.GetAtoms()):
        pos = conformer.GetAtomPosition(i)
        positions.append([pos.x, pos.y, pos.z])
        atoms.append(atom.GetSymbol())
    return np.asarray(positions), atoms

def rmsd(mol1, mol2):
    """
    Compute the RMSD between two molecules.
    """
    if isinstance(mol1, str):
        mol1 = read_molecule(mol1, remove_hs=False)
    if isinstance(mol2, str):
        mol2 = read_molecule(mol2, remove_hs=False)
    mol1,mol2 = Chem.RemoveHs(mol1), Chem.RemoveHs(mol2)
    try:
        rmsd = rdMolAlign.CalcRMS(mol1, mol2)
    except RuntimeError:
        rmsd = 100
    return rmsd

if __name__ == "__main__":
    base = "/home/haotiant/Projects/CMU/MolecularDocking"
    mol1 = f'{base}/datasets/PDBBind_processed/1a0q/1a0q_ligand.sdf'
    mol2 = f'{base}/datasets/PDBBind_processed/1a0q/1a0q_ligand_rdkit.sdf'
    print(rmsd(mol1, mol2))