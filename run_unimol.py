import os
import sys
import numpy as np
import pandas as pd
import biopandas
import lmdb
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.cluster import KMeans
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdMolAlign import AlignMolConformers
from unimol.utils.docking_utils import docking_data_pre, ensemble_iterations
from tqdm import tqdm
import pickle
import re
import json
import copy
import logging
logger = logging.getLogger(__name__)

CASF_PATH = "CASF-2016"
main_atoms = ["N", "CA", "C", "O", "H"]

def extract_pocket(pdb,ligand_coordinates,distance = 10):
    """Extract the atoms of the pocket from the protein-ligand complex by
    searching for the atoms within a certain distance from the ligand, and then extract
    the residues of found atoms.
    Args:
        pdb: str, the path to the PDB file or PDB ID.
        ligand_coordinates: 3xN array, the 3D coordinates of the ligand.
        distance: float, in the unit of Ångstrom, the distance threshold to define the pocket.
                Refernece values: 8.5 Å for water bridges, ~3.5 Å for hydrogen bonds, ~ 4.5 Å for hydrophobic interactions, ~4 Å for van der Waals interactions.

    """
    if len(pdb) == 4:
        pdb_id = pdb
        try:
            ppdb = PandasPdb().fetch_pdb(pdb_id)
        except AttributeError:
            logger.error(f"Failed to fetch the PDB ID {pdb_id}")

    elif os.path.isfile(pdb):
        pdb = os.path.abspath(pdb)
        ppdb = PandasPdb().read_pdb(pdb)
    else:
        raise ValueError("Invalid input for pdb: {}".format(pdb))
    # Extract the atoms within the distance threshold from the ligand
    coordinates = ppdb.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].values
    # Calculate the distance between each atom and the ligand atoms
    mask_protein, mask_ligand = within(coordinates, ligand_coordinates, distance)
    # combine chain and residue number to get unique residue
    ppdb.df['ATOM']["residue_chain"] = ppdb.df['ATOM']["chain_id"] + ppdb.df['ATOM']["residue_number"].astype(str)
    selected_atoms = ppdb.df['ATOM'][mask_protein]
    # Extract the residues of the selected atoms and chains
    residues = selected_atoms['residue_chain'].unique()
    # Extract the atoms of the pocket
    pocket_atoms = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_chain'].isin(residues)]
    return pocket_atoms, residues

def within(src_coords, tgt_coords, distance):
    """Check if the source coordinates are within the distance threshold from the target coordinates.
    Args:
        src_coords: Nx3 array, the source coordinates.
        tgt_coords: Mx3 array, the target coordinates.
        distance: float, the distance threshold.
    Returns:
        mask: NxM boolean array, the mask of whether the source coordinates are within the distance threshold from the target coordinates.
    """
    src_coords = src_coords[:, np.newaxis, :]
    tgt_coords = tgt_coords[np.newaxis, :, :]
    distances = np.linalg.norm(src_coords - tgt_coords, axis=-1)
    mask = distances <= distance
    mask_src = np.any(mask, axis=1)
    mask_tgt = np.any(mask, axis=0)
    return mask_src, mask_tgt

def sdf_to_dataframe(sdf_file):
    """
    Reads an SDF file and returns a pandas DataFrame with atom details.

    :param sdf_file: str, path to the SDF file
    :return: pandas DataFrame with columns 'index', 'atom_name', 'x_coord', 'y_coord', 'z_coord'
    """
    # Create an SDMolSupplier object to read the SDF file
    try:
        supplier = Chem.SDMolSupplier(sdf_file)
    except OSError:
        logger.error(f"Failed to read the SDF file {sdf_file}")
        return None, None

    # Initialize a list to hold atom data
    atom_data = []
    smiles = None
    # Iterate over all molecules in the SDF file
    for mol in supplier:
        if mol is not None:  # Check if the molecule is successfully read
            mol = Chem.RemoveHs(mol)
            try:
                smiles = Chem.MolToSmiles(mol)
            except:
                pass
            mol = AllChem.AddHs(mol, addCoords=True) #This would add potential missing hydroten in the structure
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx()
                element = atom.GetSymbol()
                pos = mol.GetConformer().GetAtomPosition(atom_idx)
                # get the SMILES
                # Add atom data to the list
                atom_data.append([atom_idx, 
                                  element, 
                                  np.float32(pos.x), 
                                  np.float32(pos.y), 
                                  np.float32(pos.z)])
            break # only read the first molecule

    # Create a DataFrame
    columns = ['index', 'atom_name', 'x_coord', 'y_coord', 'z_coord']
    df_atoms = pd.DataFrame(atom_data, columns=columns)
    print(smiles)
    return df_atoms, smiles

def load_with_ref(pdb_f, smiles, ref_ligand):
    ligand_df,smiles = sdf_to_dataframe(ref_ligand)
    pocket_atoms,pocket_residues = extract_pocket(pdb_f, ligand_df[['x_coord', 'y_coord', 'z_coord']].values)
    pmol = PandasPdb().read_pdb(pdb_f)
    return pmol, pocket_residues

def load_from_CASF(pdb_id):
    try:
        pdb_path = os.path.join(CASF_PATH, "casf2016", pdb_id + "_protein.pdb")
        pmol = PandasPdb().read_pdb(pdb_path)
        pocket_residues = json.load(
            open(os.path.join(CASF_PATH, "casf2016.pocket.json"))
        )[pdb_id]
        return pmol, pocket_residues
    except:
        print("Currently not support parsing pdb and pocket info from local files.")


def normalize_atoms(atom):
    return re.sub("\d+", "", atom)


def single_conf_gen(tgt_mol, num_confs=1000, seed=42, removeHs=True):
    mol = copy.deepcopy(tgt_mol)
    mol = Chem.AddHs(mol)
    allconformers = AllChem.EmbedMultipleConfs(
        mol, numConfs=num_confs, randomSeed=seed, clearConfs=True
    )
    sz = len(allconformers)
    for i in range(sz):
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=i)
        except:
            continue
    if removeHs:
        mol = Chem.RemoveHs(mol)
    return mol


def clustering_coords(mol, M=1000, N=100, seed=42, removeHs=True):
    rdkit_coords_list = []
    rdkit_mol = single_conf_gen(mol, num_confs=M, seed=seed, removeHs=removeHs)
    noHsIds = [
        rdkit_mol.GetAtoms()[i].GetIdx()
        for i in range(len(rdkit_mol.GetAtoms()))
        if rdkit_mol.GetAtoms()[i].GetAtomicNum() != 1
    ]
    ### exclude hydrogens for aligning
    AlignMolConformers(rdkit_mol, atomIds=noHsIds)
    sz = len(rdkit_mol.GetConformers())
    for i in range(sz):
        _coords = rdkit_mol.GetConformers()[i].GetPositions().astype(np.float32)
        rdkit_coords_list.append(_coords)

    ### exclude hydrogens for clustering
    rdkit_coords_flatten = np.array(rdkit_coords_list)[:, noHsIds].reshape(sz, -1)
    ids = (
        KMeans(n_clusters=N, random_state=seed)
        .fit_predict(rdkit_coords_flatten)
        .tolist()
    )
    coords_list = [rdkit_coords_list[ids.index(i)] for i in range(N)]
    return coords_list


def parser(pdb_id, smiles, seed=42, ref_ligand = None):
    if ref_ligand is None:
        pmol, pocket_residues = load_from_CASF(pdb_id)
    else:

    pname = pdb_id
    pro_atom = pmol.df["ATOM"]
    pro_hetatm = pmol.df["HETATM"]

    pro_atom["ID"] = pro_atom["chain_id"].astype(str) + pro_atom[
        "residue_number"
    ].astype(str)
    pro_hetatm["ID"] = pro_hetatm["chain_id"].astype(str) + pro_hetatm[
        "residue_number"
    ].astype(str)

    pocket = pd.concat(
        [
            pro_atom[pro_atom["ID"].isin(pocket_residues)],
            pro_hetatm[pro_hetatm["ID"].isin(pocket_residues)],
        ],
        axis=0,
        ignore_index=True,
    )

    pocket["normalize_atom"] = pocket["atom_name"].map(normalize_atoms)
    pocket = pocket[pocket["normalize_atom"] != ""]
    patoms = pocket["atom_name"].apply(normalize_atoms).values.tolist()
    pcoords = [pocket[["x_coord", "y_coord", "z_coord"]].values]
    side = [0 if a in main_atoms else 1 for a in patoms]
    residues = (
        pocket["chain_id"].astype(str) + pocket["residue_number"].astype(str)
    ).values.tolist()

    # generate ligand conformation
    M, N = 100, 10
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    latoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    holo_coordinates = [mol.GetConformer().GetPositions().astype(np.float32)]
    holo_mol = mol
    coordinate_list = clustering_coords(mol, M=M, N=N, seed=seed, removeHs=False)
    mol_list = [mol] * N

    return pickle.dumps(
        {
            "atoms": latoms,
            "coordinates": coordinate_list,
            "mol_list": mol_list,
            "pocket_atoms": patoms,
            "pocket_coordinates": pcoords,
            "side": side,
            "residue": residues,
            "holo_coordinates": holo_coordinates,
            "holo_mol": holo_mol,
            "holo_pocket_coordinates": pcoords,
            "smi": smiles,
            "pocket": pname,
        },
        protocol=-1,
    )


def write_lmdb(pdb_id, smiles_list, seed=42, result_dir="./results"):
    os.makedirs(result_dir, exist_ok=True)
    outputfilename = os.path.join(result_dir, pdb_id + ".lmdb")
    try:
        os.remove(outputfilename)
    except:
        pass
    env_new = lmdb.open(
        outputfilename,
        subdir=False,
        readonly=False,
        lock=False,
        readahead=False,
        meminit=False,
        max_readers=1,
        map_size=int(10e9),
    )
    for i, smiles in enumerate(smiles_list):
        inner_output = parser(pdb_id, smiles, seed=seed)
        txn_write = env_new.begin(write=True)
        txn_write.put(f"{i}".encode("ascii"), inner_output)
    txn_write.commit()
    env_new.close()


# @title Run Uni-Mol Binding Pose Prediction

# @markdown Currently this scripts only support CASF-2016 dataset with given pockets residues.

# @markdown You can input multiple SMILES, split by ','.

# @markdown If SMILES is not given, the default one in the complex will be used.

pdb_id = "4ty7"  # @param {type:"string"}
pdb_id = pdb_id.lower()
casf_collect = os.listdir(os.path.join(CASF_PATH, "casf2016"))
casf_collect = list(set([item[:4] for item in casf_collect]))
if pdb_id not in casf_collect:
  warning_str = "{} is not int CASF-2016 dataset, Please select from \n".format(pdb_id)
  for i in range(15):
    warning_str += "{}\n".format(','.join(casf_collect[20*i:20*(i+1)]))
  raise Exception(warning_str)
supp = Chem.SDMolSupplier(os.path.join(CASF_PATH, "casf2016", pdb_id + "_ligand.sdf"))
mol = [mol for mol in supp if mol][0]
ori_smiles = Chem.MolToSmiles(mol)
smiles = ""  # @param {type:"string"}
seed = 42  # @param {type:"number"}
data_path = "./CASF-2016"
results_path = "./results/"
weight_path = "/content/binding_pose_220908.pt"
batch_size = 8
dist_threshold = 8.0
recycling = 3
if smiles.split(",") == 0 or smiles == "":
    print("No other smiles inputs")
    smiles_list = [ori_smiles]
else:
    print("Docking with smiles: {}".format(smiles))
    smiles_list = smiles.split(",")

write_lmdb(pdb_id, smiles_list, seed=seed, result_dir=data_path)

!python ./Uni-Mol/unimol/unimol/infer.py --user-dir ./Uni-Mol/unimol/unimol $data_path --valid-subset $pdb_id \
       --results-path $results_path \
       --num-workers 8 --ddp-backend=c10d --batch-size $batch_size \
       --task docking_pose --loss docking_pose --arch docking_pose \
       --path $weight_path \
       --fp16 --fp16-init-scale 4 --fp16-scale-window 256 \
       --dist-threshold $dist_threshold --recycling $recycling \
       --log-interval 50 --log-format simple

def generate_docking_input(
    predict_file, reference_file, tta_times=10, output_dir="./results"
):
    (
        mol_list,
        smi_list,
        pocket_list,
        pocket_coords_list,
        distance_predict_list,
        holo_distance_predict_list,
        holo_coords_list,
        holo_center_coords_list,
    ) = docking_data_pre(reference_file, predict_file)
    iter = ensemble_iterations(
        mol_list,
        smi_list,
        pocket_list,
        pocket_coords_list,
        distance_predict_list,
        holo_distance_predict_list,
        holo_coords_list,
        holo_center_coords_list,
        tta_times=tta_times,
    )
    for i, content in enumerate(iter):
        pocket = content[3]
        output_name = os.path.join(output_dir, "{}.{}.pkl".format(pocket, i))
        try:
            os.remove(output_name)
        except:
            pass
        pd.to_pickle(content, output_name)


predict_file = os.path.join(results_path, "content_" + pdb_id + ".out.pkl")
reference_file = os.path.join(data_path, pdb_id + ".lmdb")
generate_docking_input(
    predict_file, reference_file, tta_times=10, output_dir=results_path
)
for i, smiles in enumerate(smiles_list):
    print("Docking {}".format(smiles))
    input_path = os.path.join(results_path, "{}.{}.pkl".format(pdb_id, i))
    ligand_path = os.path.join(results_path, "docking.{}.{}.sdf".format(pdb_id, i))
    cmd = "python ./Uni-Mol/unimol/unimol/utils/coordinate_model.py --input {} --output-ligand {}".format(
        input_path, ligand_path
    )
    os.system(cmd)