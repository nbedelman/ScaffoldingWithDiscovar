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
module load bedtools
export PATH=$PATH:/n/mallet_lab/edelman/software/kentUtils/bin/

round=$1
chromosomeResultsDir=$2
outputBase=$3
overviewBedFile=$4

#find scaffolds where a large part was deleted
cat $chromosomeResultsDir/*coordinates.txt > $outputBase\_allCoordinates.txt
coordinatesToTSV.py $outputBase\_allCoordinates.txt $outputBase\_allCoordinates.tsv
findSuspiciousScaffolds.py data/refGenomeContigs.size $outputBase\_allCoordinates.tsv $outputBase\_highlyAlteredScaffolds.csv

#combine results into a single chromosome-level fasta file
makeFullChromFasta.sh $round $chromosomeResultsDir $outputBase

#split into contigs and get stats
makeGenomeMap.py $outputBase\_chroms.fa 100 name
mv $outputBase\_chroms.fa.bed $outputBase\_allContigs.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_allContigs.bed -fo $outputBase\_allContigs.fa
faSize $outputBase\_allContigs.fa -detailed > $outputBase\_allContigs.size
grep -v Ns  $outputBase\_allContigs.bed > $outputBase\_noNs_contigs.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_noNs_contigs.bed -fo $outputBase\_noNs_contigs.fa
genomeStats.py $outputBase\_noNs_contigs.bed bed $outputBase\_noNs_contigs.stats

#split into scaffolds and get stats
makeGenomeMap.py $outputBase\_chroms.fa 1000 name exact
mv $outputBase\_chroms.fa.bed $outputBase\_allScaffolds.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_allScaffolds.bed -fo $outputBase\_allScaffolds.fa
faSize $outputBase\_allScaffolds.fa -detailed > $outputBase\_allScaffolds.size
grep -v Ns $outputBase\_allScaffolds.bed > $outputBase\_noNs_scaffolds.bed
bedtools getfasta -s -name -fi $outputBase\_chroms.fa -bed $outputBase\_noNs_scaffolds.bed -fo $outputBase\_noNs_scaffolds.fa
genomeStats.py $outputBase\_noNs_scaffolds.bed bed $outputBase\_noNs_scaffolds.stats

#get mapping quality analysis
mappingQualityAnalysis.R $overviewBedFile

#move results into the results folder.
mv $outputBase* overallResults
mv *.pdf overallResults
