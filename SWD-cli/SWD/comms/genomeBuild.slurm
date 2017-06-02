#!/bin/bash

#genomeBuild.slurm
#runs the lastdb command

#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-02:00 # Runtime
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=5000 # Memory per cpu in MB
#SBATCH -o lastdb.out # Standard out goes to this file
#SBATCH -e lastdb.err # Standard err goes to this file

options=$1
dbName=$2
refGenome=$3

lastdb $options $dbName $refGenome
mkdir -p lastdb
mv $dbName.* lastdb/
