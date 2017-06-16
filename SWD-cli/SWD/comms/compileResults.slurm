#!/bin/bash

#compileResults.slurm

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=1000
#SBATCH -t 0-00:30
#SBATCH -e compileResults.err
#SBATCH -o compileResults.out

module load R

round=$1
chromosomeResultsDir=$2
outputBase=$3
overviewBedFile=$4

cat $chromosomeResultsDir/*coordinates.txt > $outputBase\_allCoordinates.txt
coordinatesToTSV.py $outputBase\_allCoordinates.txt $outputBase\_allCoordinates.tsv
makeFullChromFasta.sh $round $chromosomeResultsDir $outputBase
makeGenomeMap.py $outputBase\_chroms.fa 100 name
mv $outputBase\_chroms.fa.bed $outputBase\_chroms.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_chroms.bed -fo $outputBase\_allContigs.fa
grep -v Ns  $outputBase\_chroms.bed > $outputBase\_noNs_contigs.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_noNs_contigs.bed -fo $outputBase\_noNs_contigs.fa
genomeStats.py $outputBase\_noNs_contigs.bed bed $outputBase\_noNs_contigs.stats
mappingQualityAnalysis.R $overviewBedFile
mv $outputBase* overallResults
mv *.pdf overallResults
