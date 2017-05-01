#!/bin/bash
#
# genomeAlign_last.slurm from tophat.sh example
# This script aligns the new assembly to the reference genome, and creates
# custom bed files for the next steps of the pipeline. It was created to run
# on a an assembly that has been split into several parts, and it runs the
# lastal alignment and mafToBed script on each part separately.

#SBATCH -J lastal         # A single job name for the array
#SBATCH -p serial_requeue # Partition
#SBATCH -n 8             # eight cores
#SBATCH -t 2-00:00         # Running time of 2 days
#SBATCH --mem 5000        # Memory request of 2 GB
#SBATCH -o outs/lastal%A_%a.out # Standard output
#SBATCH -e errs/lastal%A_%a.err # Standard error

reference=$1
speciesLoc=$2
species=$3
lastalConfig=$4
discoData=$5

mkdir -p $species
mkdir -p $species/mafs
mkdir -p $species/beds
mkdir -p $species/sams

refLoc=$reference
query=$speciesLoc/parts/*.part-${SLURM_ARRAY_TASK_ID}.fa
alignOut=$species/mafs/$species.part-${SLURM_ARRAY_TASK_ID}_alns.maf
samOut=$species/sams/$species.part-${SLURM_ARRAY_TASK_ID}_alns.sam


#find local alignments between query sequences and reference sequences
SWD lastAlign $refLoc $query $lastalConfig | maf-sort.sh > $alignOut
#maf-convert.py sam $alignOut > $samOut

##Make a folder in which to put the bed files for each contig; then make the files
mkdir $species/beds/${SLURM_ARRAY_TASK_ID}_tmp
cd $species/beds/${SLURM_ARRAY_TASK_ID}_tmp

SWD mafToBed ../../../$alignOut ../../../$discoData
#mafToBed.py ../../../$alignOut $discoData

#put a copy of each bed file in the general bed folder, and concatenate into the allContigs file.
#Also put the overviews into the allOverviews file. Then get rid of the temporary folder.
for i in *.bed
do cp $i ../$i
cat $i >> ${SLURM_ARRAY_TASK_ID}.bed
tail -n +2 $i >> ../$species\_allContigs.bed
head -n 1 $i >> ../$species\_allOverviews.bed
done

cd ..
rm -Rf ${SLURM_ARRAY_TASK_ID}_tmp
