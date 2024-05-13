#!/bin/sh

#some directories that will never be in version control
touch .gitignore
echo "conda-env" > .gitignore
echo "workflow/.snakemake" >> .gitignore
echo "workflow/slurm-*.out" >> .gitignore
git add .gitignore

#some boilerplate for the workflow
mkdir -p workflow/envs
mkdir -p workflow/scripts
mkdir -p workflow/reports

#make the Snakemake header
echo "###############################################################################
#Snakemake
#About:
#
#
#
#@author: charles.hefer@agresearch.co.nz
###############################################################################" > workflow/Snakefile


#some config related files
mkdir -p ./config
echo "#!/bin/sh
conda create -p ./conda-env -c conda-forge -c bioconda -y snakemake==7.32.4 mamba python apptainer --solver libmamba" > config/setup.sh
cp ~/Templates/config.yaml config/config.yaml

#the status checking file for slurm
cp ~/Templates/status.py config/status.py

#The slurm config
mkdir -p ./slurm
cp ~/Templates/slurm_config.yaml slurm/config.yaml

mkdir -p ./results
mkdir -p ./resources

