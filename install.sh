# this script install the specified package with input arguments
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi
#package including diffbindfr, diffdock, unimol, dtmol, smina and gnina

#install diffbindfr
function install_diffbindfr {
    echo "Installing diffbindfr..."
    git clone https://github.com/HBioquant/DiffBindFR.git
    cd DiffBindFR
    #check if cuda version is 11.7
    if [ $(nvcc --version | grep "release 11.7" | wc -l) -ne 1 ]; then
        echo "DiffBindFR requires CUDA 11.7. Please switch CUDA to 11.7 and try again."
        exit 1
    fi
    #pip uninstall -y torch-scatter torch-cluster torch-sparse torch-spline-conv torchmetrics torch-geometric
    conda env create -f env.yaml
    conda activate diffbindfr
    pip install -e .
    conda deactivate
    cd ..
    echo "DiffBindFR has been installed."
}

#install diffdock
function install_diffdock {
    echo "Installing diffdock..."
    git clone https://github.com/gcorso/DiffDock.git
    cd DiffDock
    git checkout 5238b18
    if [ $(nvcc --version | grep "release 11.7" | wc -l) -ne 1 ]; then
        echo "DiffBindFR requires CUDA 11.7. Please switch CUDA to 11.7 and try again."
        exit 1
    fi
    conda env create --file environment.yml
    conda activate diffdock
    pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
    pip install torch_scatter torch_sparse torch_cluster -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
    pip install prody rdkit networkx fair-esm[esmfold]==2.0.0 e3nn==0.5.1 
    cd data
    wget https://zenodo.org/records/10656052/files/BindingMOAD_2020_processed.tar?download=1
    wget https://zenodo.org/records/10656052/files/DockGen.tar?download=1

}

#install unimol
function install_unimol {
    echo "Installing unimol..."
    #git clone
    #TODO
}

#install dtmol
function install_dtmol {
    echo "Installing dtmol..."
    #git clone
    #TODO
}

#install smina
function install_smina {
    echo "Installing smina..."
    mkdir -p smina
    cd smina
    wget https://sourceforge.net/projects/smina/files/smina.static/download -O smina
    chmod +x smina
    #TODO
}

#install gnina
function install_gnina {
    echo "Installing gnina..."
    mkdir -p gnina
    cd gnina
    wget https://github.com/gnina/gnina/releases/download/v1.1/gnina
    chmod +x gnina

    # or install from source (require sudo permission to successfully install)
    # BASE=$(pwd)
    # mkdir -p run
    # git clone https://github.com/openbabel/openbabel.git
    # cd openbabel
    # mkdir build
    # cd build
    # cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DCMAKE_INSTALL_PREFIX=${BASE}/run/ ..
    # make -j4
    # make install

    # cd ${BASE}
    # wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-Linux-x86_64.tar.gz
    # tar -xvf cmake-3.25.2-Linux-x86_64.tar.gz
    # rm cmake-3.25.2-Linux-x86_64.tar.gz
    # if [ $(nvcc --version | grep "release 12" | wc -l) -ne 1 ]; then
    #     echo "gnina requires CUDA > 12.0. Please switch CUDA to 12.0 or higher and try again."
    #     exit 1
    # fi
    # git clone https://github.com/gnina/gnina.git
    # cd gnina
    # mkdir build 
    # cd build
    # cmake -DCMAKE_PREFIX_PATH=${BASE}/run/ -DOpenBabel3_DIR=${BASE}/run/ ..
    # make -j4
    # make install

}

#install
for arg in $@; do
    case $arg in
        "diffbindfr")
            install_diffbindfr
            ;;
        "diffdock")
            install_diffdock
            ;;
        "unimol")
            install_unimol
            ;;
        "dtmol")
            install_dtmol
            ;;
        "smina")
            install_smina
            ;;
        "gnina")
            install_gnina
            ;;
        *)
            echo "Invalid package name: $arg"
            ;;
    esac
done
