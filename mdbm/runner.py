#This script provide general method to run pdbbind_test dataset benchamrk given a method
import os
import numpy as np
from typing import Callable
from rdkit import Chem
from mdbm.utils.mol import rmsd,read_molecule
from time import time

class Timer(object):
    def __init__(self):
        self.event_list =  None
        self.curr_event = -1
        self.curr_time = time()
        self.records = {}
        self.timing = False

    def reset_for(self, events):
        """Reset the timer and set the count_list
        events: A list of events need to record the time
        """
        self.event_list = events
        self.curr_event = -1
        self.timing = False

    def start_next(self):
        self.curr_event += 1
        assert not self.timing, "You need to stop the current timing event first"
        self.timing = True
        self.curr_time = time()

    def end_current(self):
        self.records[self.event_list[self.curr_event]] = time() - self.curr_time
        assert self.timing, "You need to start the timing event first"
        self.timing = False

    def average_time(self):
        return np.mean(list(self.records.values()))
    

class Runner(object):
    def __init__(self, df, save_dir, mode = "nohup"):
        """
        df: str, path to the directory of the pdbbind_test dataset
        save_dir: str, path to the directory to save the results
        mode: str, the mode to run, can be 'nohup' mode and 'debug' mode, will first didn't stop running will fail and keep running.
        """
        self.df = df
        self.save_dir = save_dir
        self.pdb_ids = list(os.listdir(df))
        self.mode = mode
        self.timer = Timer()

    def run(self,method:Callable):
        """
        method: Callable, the method to run on the dataset, 
            requires 3 or 4 arguments: receptor_f, ligand_f, ref_l, out 
            returns a molecule object
        """
        raise NotImplementedError
    
    def _run(self,method:Callable, protein_postfix, ref_ligand_postfix, ligand_postfix):
        rmsd_vals = []
        scores_vals = []
        pdb_ids = []
        #empty the report file
        report_f = os.path.join(self.save_dir, 'rmsd_report.txt')
        with open(report_f, 'w') as f:
            f.write("PDB ID\tRMSD(conformers)\tElapsed Time\n")
        self.timer.reset_for(self.pdb_ids)
        for pdb_id in self.pdb_ids:
            self.timer.start_next()
            receptor_f = os.path.join(self.df, pdb_id, f"{pdb_id}{protein_postfix}")
            ligand_f = os.path.join(self.df, pdb_id, f"{pdb_id}{ref_ligand_postfix}")
            ref_l = os.path.join(self.df, pdb_id, f"{pdb_id}{ligand_postfix}")
            out = os.path.join(self.save_dir, f"{pdb_id}_result.sdf")
            try:
                mol,scores = method(receptor_f, ligand_f, ref_l, out)
                curr_rmsd = self.get_rmsd(mol, read_molecule(ref_l))
            except Exception as e:
                if self.mode == "nohup":
                    print(f"Failed to run {pdb_id}: {e}")
                    self.timer.end_current()
                    continue
                else:
                    raise e
            self.timer.end_current()
            #rank the rmsd values by the scores
            curr_rmsd = [x for _,x in sorted(zip(scores, curr_rmsd))]
            rmsd_vals.append(curr_rmsd)
            scores = sorted(scores)
            scores_vals.append(scores)
            pdb_ids.append(pdb_id)
            self.single_rmsd_report(pdb_id, curr_rmsd, scores)
        return self.rmsd_reports(pdb_ids, rmsd_vals, scores_vals)

    def get_rmsd(self, mol_cluster, target_mol):
        """
        Get the top n rmsd values from the mol_cluster
        Args:
            mol_cluster: list of molecule objects
            target_mol: molecule object
            scores: A list of scores for the molecules
        """
        rmsd_vals = []
        for mol in mol_cluster:
            rmsd_vals.append(rmsd(mol, target_mol))
        return sorted(rmsd_vals)
    
    def single_rmsd_report(self,pdb_id, rmsd_vals, scores):
        report_f = os.path.join(self.save_dir, 'rmsd_report.txt')
        with open(report_f, 'a') as f:
            msg = f"{pdb_id}\t{' '.join([str(x) for x in rmsd_vals])}\t{self.timer.records[pdb_id]}\n"
            f.write(msg)
            print(msg)

    def rmsd_reports(self,pdbs, rmsd_vals, scores):
        """
        Generate a well-formated report of the rmsd values
        Args:
            pdbs: list of pdb ids
            rmsd_vals: list of list of rmsd values
            scores: list of list of scores
        """
        summary_f = os.path.join(self.save_dir, 'rmsd_summary.txt')
        with open(summary_f, 'w') as f:
            #report porpotion of rmsd < 2
            f.write(f"Proportion of RMSD < 2 (Top 1): {len([r[0] for r in rmsd_vals if r[0] < 2])/len(rmsd_vals)}")
            #report porpotion of rmsd < 1
            f.write(f"Proportion of RMSD < 1 (Top 1): {len([r[0] for r in rmsd_vals if r[0] < 1])/len(rmsd_vals)}")
            #report porpotion of rmsd < 2 for top 5 scores rmsd
            f.write(f"Proportion of RMSD < 2 (Top 5): {len([min(r[:5]) for r in rmsd_vals if min(r[:5]) < 2])/len(rmsd_vals)}")
            #report propotion of rmsd <1 for top 5 scores rmsd
            f.write(f"Proportion of RMSD < 1 (Top 5): {len([min(r[:5]) for r in rmsd_vals if min(r[:5]) < 1])/len(rmsd_vals)}")
            #report the median top1 rmsd
            f.write(f"Median RMSD (Top 1): {np.median([r[0] for r in rmsd_vals])}")
            #report the median top5 rmsd
            f.write(f"Median RMSD (Top 5): {np.median([min(r[:5]) for r in rmsd_vals])}")
            #report for the running time:
            f.write(f"Average running time {self.timer.average_time()}")

        return summary_f

class PDBBindRunner(Runner):
    def run(self,method:Callable):
        """
        Run the method on the pdbbind_test dataset
        """
        protein_postfix = '_protein_processed.pdb'
        ref_ligand_postfix = '_ligand.sdf'
        ligand_postfix = '_ligand.sdf'
        return self._run(method, protein_postfix, ref_ligand_postfix, ligand_postfix)


class PBRunner(Runner):
    def run(self,method:Callable):
        """
        Run the method on the pb dataset
        """
        protein_postfix = '_protein.pdb'
        ref_ligand_postfix = '_ligand.sdf'
        ligand_postfix = '_ligand_start_conf.sdf'
        return self._run(method, protein_postfix, ref_ligand_postfix, ligand_postfix)