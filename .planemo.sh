#!/bin/bash

tmp_dir=`mktemp -d`
source activate python2.7
# planemo_bin='which planemo'
echo "GOOO"
# echo $planemo_bin
echo "BLA"
# source deactivate
planemo database_create galaxy
planemo conda_init --conda_prefix $tmp_dir/conda
export PATH=$tmp_dir/conda/bin:$PATH
conda install -y -c bioconda samtools python=2.7.13 numpy scipy matplotlib=2.0.0 nose flake8 pytables biopython pysam pybigwig intervaltree future six pandas

# source activate hicexplorer_galaxy

pip install .


# Galaxy wrapper testing
planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
# /home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
source deactivate