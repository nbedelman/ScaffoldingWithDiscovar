#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --mem=50000
#SBATCH -o repeatMask.out
#SBATCH -e repeatMask.err
#SBATCH --constraint holyib

genome=$1
outDir=$2

RepeatMasker -pa 32 -species hexapoda -xsmall -nocut -dir $outDir -q $genome
