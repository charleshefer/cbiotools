#!/bin/bash -e

#SBATCH --job-name=test       #A sensible name for this job
#SBATCH --ntasks=1              # Num of processes can be lanuched by this job
#SBATCH --cpus-per-task=1       # Num of CPUs allocated to each lanuched process
#SBATCH --mem=2G
#SBATCH --time=00:01:00         # The allowed duration of this job in hh:mm:ss format
#SBATCH --partition=compute,hugemem  # Use the compute and hugemem partition
#SBATCH --mail-user
#SBATCH --mail-type=END
#SBATCH --account=2014-white_clover_gbs

sleep 10

echo $HOSTNAME > hostname.txt



