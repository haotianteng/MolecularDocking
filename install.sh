# this script install the specified package with input arguments
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi
#package including diffbindfr, diffdock, unimol, dtmol, smina and gnina

#install diffbindfr
function install_diffbindfr {
    echo "Installing diffbindfr..."
    #git clone https://github.com/HBioquant/DiffBindFR.git
    cd DiffBindFR
    #check if cuda version is 11.7
    if [ $(nvcc --version | grep "release 11.7" | wc -l) -ne 1 ]; then
        echo "DiffBindFR requires CUDA 11.7. Please switch CUDA to 11.7 and try again."
        exit 1
    fi
    #pip uninstall -y torch-scatter torch-cluster torch-sparse torch-spline-conv torchmetrics torch-geometric
    conda env create -f diffbindfr.yml
    conda activate diffbindfr
    pip install -e .
    conda deactivate
    cd ..
    echo "DiffBindFR has been installed."
}
    