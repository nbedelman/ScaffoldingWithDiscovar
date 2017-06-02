#!/bin/bash

#reOrderScaffolds.slurm
#runs the SWD reOrderScaffolds command

#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-05:00 # Runtime of five hours
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=10000 # Memory per cpu in MB
#SBATCH -o reOrderScaffolds.out # Standard out goes to this file
#SBATCH -e reOrderScaffolds.err # Standard err goes to this file

refGenome=$1
newAssembly=$2
ungroupedChrom=$3
bedDirectory=$4
map=$5

SWD.py reOrder --refGenome $refGenome --newAssembly $newAssembly --ungroupedChrom $ungroupedChrom $bedDirectory $map

mkdir -p chromosomeResults
mv *discoOrder* chromosomeResults
