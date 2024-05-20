#Script to run smina on pb and PDBBind datasets
import argparse
import os
import sys
from rdkit import Chem
from functools import partial
from mdbm.utils.mol import read_molecule
from mdbm.runner import PDBBindRunner, PBRunner

CURR_F= os.path.dirname(os.path.abspath(__file__))

def run_once(receptor_f, ligand_f, orig_ligand_f, out_f, box_add=5, **kwargs):
    additional_args = " ".join([f"{k} {v}" for k,v in kwargs.items()])
    cmd = f'./smina/smina -r {receptor_f} -l {ligand_f} --autobox_ligand {orig_ligand_f} -o {out_f} --autobox_add {box_add} {additional_args}'
    print(cmd)
    os.system(cmd)
    #read back the output file and return the mol
    mol = Chem.SDMolSupplier(out_f)
    scores = [float(m.GetProp('minimizedAffinity')) for m in mol]


def run_pb(df, out_f, mode, box_add=5, **kwargs):
    runner = PBRunner(df, out_f, mode=mode)
    method = partial(run_once, box_add=box_add, **kwargs)
    runner.run(method)

def run_pdbbind(df, out_f, mode, box_add=5, **kwargs):
    runner = PDBBindRunner(df, out_f, mode=mode)
    method = partial(run_once, box_add=box_add, **kwargs)
    runner.run(method)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run smina')
    parser.add_argument('--dataset', type=str, default='pb', help='dataset name to run')
    parser.add_argument('--box_add', type=int, default=5, help='box add value')
    parser.add_argument('--mode', type=str, default='nohup', help='mode to run, can be nohup or debug')
    #take the rest of the arguments as smina arguments
    args,additional_args = parser.parse_known_args(sys.argv[1:])
    smina_args = {}
    for i in range(0,len(additional_args),2):
        smina_args[additional_args[i]] = additional_args[i+1]
    out_f = os.path.join(CURR_F, 'output', 'smina_out')
    os.makedirs(out_f, exist_ok=True)
    if args.dataset == 'pb':
        out_f1 = os.path.join(out_f, 'pb_astex')
        os.makedirs(out_f1, exist_ok=True)
        run_pb(os.path.join(CURR_F, 'datasets/pb/astex_diverse_set'), out_f1, mode = args.mode, box_add=args.box_add, **smina_args)
        out_f2 = os.path.join(out_f, 'pb')
        os.makedirs(out_f2, exist_ok=True)
        run_pb(os.path.join(CURR_F, 'datasets/pb/posebusters_benchmark_set'), out_f2, mode = args.mode, box_add=args.box_add, **smina_args)
    elif args.dataset == 'pdbbind':
        run_pdbbind(os.path.join(CURR_F, 'datasets/PDBBind_processed_test'), 
                    out_f, 
                    mode = args.mode,
                    box_add=args.box_add)
    elif args.dataset == 'test':
        run_pb(os.path.join(CURR_F, 'datasets/test'), out_f, mode = args.mode, use_cnn_score = args.use_cnn_score, box_add=args.box_add, **smina_args)
    else:
        print('Invalid dataset name')
        sys.exit(1)

