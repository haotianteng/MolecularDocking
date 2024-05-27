### This sciprt is used to download specify datasets from the internet
#download from https://zenodo.org/records/6408497/files/PDBBind.zip
#Use multiple threads
mkdir -p datasets
cd datasets
if [ ! -d "PDBBind_processed" ]; then
    wget https://zenodo.org/records/6408497/files/PDBBind.zip
    unzip PDBBind.zip
    rm PDBBind.zip
else
    echo "PDBBind data folder already exists"
fi

#download posebuster from to pb folder
if [ ! -d "pb" ]; then
    mkdir -p pb
    wget https://zenodo.org/records/8278563/files/posebusters_paper_data.zip -O pb/posebusters_paper_data.zip
    cd pb
    unzip posebusters_paper_data.zip
    rm posebusters_paper_data.zip
    cd ..
else
    echo "posebuster data folder already exists"
fi

#split the pdbbind dataset by ../pdbbind_timesplit.txt
if [ ! -d "PDBBind_processed_test" ]; then
    #read in the pdb_ids in ../pdbbind_timesplit.txt file
    mkdir -p PDBBind_processed_test
    while IFS= read -r line
    do
        echo "Copying $line to PDBBind_processed_test"
        #copy the folder to PDBBind_processed_test
        cp -r PDBBind_processed/$line PDBBind_processed_test/
    done < ../pdbbind_timesplit.txt

else
    echo "PDBBind_processed_test folder already exists"
fi

#Check if PDBBind dataset is downloaded correctly and has 19120 folders in PDBBind_processed
echo "Checking if PDBBind dataset is downloaded correctly..."
if [ ! -d "PDBBind_processed" ]; then
    echo "PDBBind dataset is not downloaded correctly. Please check the download link."
    exit 1
else
    if [ $(ls PDBBind_processed | wc -l) -ne 19120 ]; then
        echo "PDBBind dataset is not downloaded correctly. Please check the download link."
        exit 1
    else
        echo "PDBBind dataset has been downloaded correctly."
    fi
fi

#check if PDBBind_processed_test has 363 folders
echo "Checking if PDBBind_processed_test dataset is splited correctly..."
if [ ! -d "PDBBind_processed_test" ]; then
    echo "PDBBind_processed_test dataset is not splited correctly. Please check the download link."
    exit 1
else
    if [ $(ls PDBBind_processed_test | wc -l) -ne 363 ]; then
        echo "PDBBind_processed_test dataset is not splited correctly. Please check the download link."
        exit 1
    else
        echo "PDBBind_processed_test dataset has been splited correctly."
    fi
fi


#Check if posebuster dataset is downloaded correctly and has 428 folders in pb/posebusters_benchmark_set and 85 folders in pb/astex_diverse_set
echo "Checking if posebuster dataset is downloaded correctly..."
if [ ! -d "pb/posebusters_benchmark_set" ]; then
    echo "posebuster dataset is not downloaded correctly. Please check the download link."
    exit 1
else
    if [ $(ls pb/posebusters_benchmark_set | wc -l) -ne 428 ]; then
        echo "posebuster benchmark_set is not downloaded correctly. Please check the download link."
        exit 1
    else
        if [ $(ls pb/astex_diverse_set | wc -l) -ne 85 ]; then
            echo "posebuster diverse_set is not downloaded correctly. Please check the download link."
            exit 1
        fi
        echo "posebuster dataset has been downloaded correctly."
    fi
fi

### Create the test dataset
mkdir -p test
cp -r pb/astex_diverse_set/1G9V_RQ3 test/
cp -r pb/astex_diverse_set/1GKC_NFH test/
cp -r pb/astex_diverse_set/1GM8_SOX test/


echo "All datasets have been downloaded."

