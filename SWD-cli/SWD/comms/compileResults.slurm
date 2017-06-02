#!/bin/bash

#compileResults.slurm

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=1000
#SBATCH -t 0-00:30
#SBATCH -e compileResults.err
#SBATCH -o compileResults.out

round=$1
chromosomeResultsDir=$2
outputBase=$3
overviewBedFile=$4


makeFullChromFasta.sh $round $chromosomeResultsDir $outputBase
makeGenomeMap.py $outputBase\_chroms.fa 100 name
mv $outputBase\_chroms.fa.bed $outputBase\_chroms.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_chroms.bed -fo $outputBase\_allContigs.fa
grep -v Ns  $outputBase\_chroms.bed > $outputBase\_noNs_contigs.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_noNs_contigs.bed -fo $outputBase\_noNs_Contigs.fa
genomeStats.py $outputBase\_noNs_contigs.bed bed $outputBase\_noNs_Contigs.stats
mappingQualityAnalysis.R $overviewBedFile 
mv $outputBase* overallResults
mv *.pdf overallResults
