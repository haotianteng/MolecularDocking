#Script to run smina on pb and PDBBind datasets
import argparse
import os
import sys
import re
from rdkit import Chem
from functools import partial
from mdbm.utils.mol import read_molecule
from mdbm.runner import PDBBindRunner, PBRunner

CURR_F= os.path.dirname(os.path.abspath(__file__))

def run(out_f):
    # This will run Diffdock on the protein ligand list files
    cmd = f"python DiffDock/inference.py --config DiffDock/default_inference_args.yaml \
        --protein_ligand_csv {out_f}/protein_ligand_list.csv --out_dir {out_f}\
        --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise"
    os.system(cmd)

def run_once(receptor_f, ligand_f, orig_ligand_f, out_f, box_add=5, **kwargs):
    #get PDB id fromt he receptor file
    #match pattern PDB_ID_result.sdf\
    pdb_id = os.path.basename(out_f)
    pdb_id = re.match(r'(\w+)_result.sdf',pdb_id)[1]
    out_f = os.path.join(os.path.dirname(out_f),pdb_id)
    mols,scores = [],[]
    for f in os.listdir(out_f):
        if ('rank' in f) and ('confidence' in f):
            mol = Chem.SDMolSupplier(os.path.join(out_f,f))
            #match the pattern rank%d_confidence%f.sdf
            print(f)
            score = float(re.match(r'rank(\d+)_confidence(-?\d+\.\d+).sdf',f)[2])
            mols.append(mol[0])
            scores.append(score)
    return mols,scores


def prepare_run_file(df,out_f,dataset = "pb"):
    protein_postfix = "_protein.pdb" if dataset == "pb" else "_protein_processed.pdb"
    ligand_postfix = "_ligand_start_conf.sdf" if dataset == "pb" else "_ligand_rdkit.sdf"
    with open(os.path.join(out_f,'protein_ligand_list.csv'),'w+') as f:
        f.write('complex_name,protein_path,ligand_description,protein_sequence\n')
        for pdb_id in os.listdir(df):
            f.write(f'{pdb_id},{os.path.join(df,pdb_id,f"{pdb_id}{protein_postfix}")},{os.path.join(df,pdb_id,f"{pdb_id}{ligand_postfix}")},\n')

def run_pb(df, out_f, mode, box_add=5, **kwargs):
    prepare_run_file(df,out_f,dataset = "pb")
    run(out_f)
    runner = PBRunner(df, out_f, mode=mode)
    method = partial(run_once, box_add=box_add, **kwargs)
    runner.run(method)

def run_pdbbind(df, out_f, mode, box_add=5, **kwargs):
    prepare_run_file(df,out_f,dataset = "pdbbind")
    run(out_f)
    runner = PDBBindRunner(df, out_f, mode=mode)
    method = partial(run_once, box_add=box_add, **kwargs)
    runner.run(method)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run smina')
    parser.add_argument('--dataset', type=str, default='test', help='dataset name to run')
    parser.add_argument('--output', type=str, default='diffdock_out', help='output directory')
    parser.add_argument('--mode', type=str, default='nohup', help='mode to run, can be nohup or debug')
    #take the rest of the arguments as smina arguments
    args,additional_args = parser.parse_known_args(sys.argv[1:])
    docking_args = {}
    for i in range(0,len(additional_args),2):
        docking_args[additional_args[i]] = additional_args[i+1]
    out_f = os.path.join(CURR_F, 'output', args.output)
    os.makedirs(out_f, exist_ok=True)
    if args.dataset == 'pb':
        out_f1 = os.path.join(out_f, 'pb_astex')
        os.makedirs(out_f1, exist_ok=True)
        run_pb(os.path.join(CURR_F, 'datasets/pb/astex_diverse_set'), out_f1, mode = args.mode, **docking_args)
        out_f2 = os.path.join(out_f, 'pb')
        os.makedirs(out_f2, exist_ok=True)
        run_pb(os.path.join(CURR_F, 'datasets/pb/posebusters_benchmark_set'), out_f2, mode = args.mode, **docking_args)
    elif args.dataset == 'pdbbind':
        out_f = os.path.join(out_f, 'pdbbind')
        os.makedirs(out_f,exist_ok = True)
        run_pdbbind(os.path.join(CURR_F, 'datasets/PDBBind_processed_test'), 
                    out_f, 
                    mode = args.mode,
        )
    elif args.dataset == 'test':
        out_f = os.path.join(out_f,'test')
        os.makedirs(out_f,exist_ok = True)
        run_pb(os.path.join(CURR_F, 'datasets/test'), out_f, mode = args.mode, **docking_args)
    else:
        print('Invalid dataset name')
        sys.exit(1)


