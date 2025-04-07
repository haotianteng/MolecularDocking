# Molecular Docking Benchmark
Repository holding the code for molecular docking benchmark

---
## Table of contents

- [Table of contents](#table-of-contents)
- [Install](#install)
- [Dataset Download](#dataset)
- [Run Benchmark](#run)
  
## Install
For a painless installation, recommend to create a vritual environment for each docking tool. For example install for smina
```bash
conda create --name smina
conda activate smina
bash install.sh smina
```
Installing DiffbindFR and Diffdock using the script will automatically create a environment for you.

## Dataset
To download the datasets (PDBBind, PoseBusters), run the download script
```bash
bash download_datasets.sh
```
This will download the PDBBind and PoseBusters to the repository.

## Run
To run evaluation for model for dataset, use the corresponding script, for example, to get the benchmark result of gnina on PDBBind
```bash
conda activate GNINA
python run_gninia.py --dataset pdbbind
```
Dataset can be `pdbbind` (PDBBind) or `pb` (PoseBuster)
