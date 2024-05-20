#Script to run gnina on pb and PDBBind datasets
import argparse
import os
import sys
from rdkit import Chem
from functools import partial
from mdbm.utils.mol import read_molecule
from mdbm.runner import PDBBindRunner, PBRunner

CURR_F= os.path.dirname(os.path.abspath(__file__))

def run_once(receptor_f, ligand_f, orig_ligand_f, out_f, box_add=5, use_cnn_score = True):
    cmd = f'./gnina/gnina -r {receptor_f} -l {ligand_f} --autobox_ligand {orig_ligand_f} -o {out_f} --autobox_add {box_add}'
    print(cmd)
    os.system(cmd)
    #read back the output file and return the mol
    mol = Chem.SDMolSupplier(out_f)
    if use_cnn_score:
        scores = [-float(m.GetProp('CNNscore')) for m in mol]
    else:
        scores = [float(m.GetProp('minimizedAffinity')) for m in mol]
    return mol, scores

def run_pb(df, out_f, box_add=5):
    runner = PBRunner(df, out_f)
    method = partial(run_once, box_add=box_add)
    runner.run(method)

def run_pdbbind(df, out_f, box_add=5):
    runner = PDBBindRunner(df, out_f)
    method = partial(run_once, box_add=box_add)
    runner.run(method)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run gnina')
    parser.add_argument('--dataset', type=str, default='pb', help='dataset name to run')
    parser.add_argument('--box_add', type=int, default=5, help='box add value')
    args = parser.parse_args(sys.argv[1:])
    out_f = os.path.join(CURR_F, 'output', 'gnina_out')
    os.makedirs(out_f, exist_ok=True)
    if args.dataset == 'pb':
        out_f1 = os.path.join(out_f, 'pb_astex')
        os.makedirs(out_f1, exist_ok=True)
        run_pb(os.path.join(CURR_F, 'datasets/pb/astex_diverse_set'), out_f1, box_add=args.box_add)
        out_f2 = os.path.join(out_f, 'pb')
        os.makedirs(out_f2, exist_ok=True)
        run_pb(os.path.join(CURR_F, 'datasets/pb/posebusters_benchmark_set'), out_f2, box_add=args.box_add)
    elif args.dataset == 'pdbbind':
        run_pdbbind(os.path.join(CURR_F, 'datasets/PDBBind_processed_test'), 
                    out_f, 
                    box_add=args.box_add)
    elif args.dataset == 'test':
        run_pb(os.path.join(CURR_F, 'datasets/test'), out_f, box_add=args.box_add)
    else:
        print('Invalid dataset name')
        sys.exit(1)

