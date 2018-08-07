#!/bin/bash

#overviewToIntersects.slurm
#runs the SWD overviewToIntersects command

#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 00:30:00 # Runtime
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=500 # Memory per cpu in MB
#SBATCH -o overviewToIntersects.out # Standard out goes to this file
#SBATCH -e overviewToIntersects.err # Standard err goes to this file

overview=$1
bedDirectory=$2

SWD.py overviewToIntersects $overview $bedDirectory 
