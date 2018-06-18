#!/bin/bash

#compileResults.slurm

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=1000
#SBATCH -t 0-00:30
#SBATCH -e compileResults.err
#SBATCH -o compileResults.out

module load Anaconda
module load R
module load bedtools
module load RepeatMasker
export PATH=$PATH:/n/mallet_lab/edelman/software/kentUtils/bin/
export PATH=$PATH:/n/mallet_lab/edelman/software/PBSuite_15.8.24/pbsuite/utils/
source activate my_root

round=$1
chromosomeResultsDir=$2
outputBase=$3
overviewBedFile=$4

#move all chromosomal results to chromosomeResultsDir
mv *coordinates.txt *report.csv *_overwrittenScaffolds.txt *_envelopers.txt $chromosomeResultsDir

#find scaffolds where a large part was deleted
cat $chromosomeResultsDir/*coordinates.txt > $outputBase\_allCoordinates.txt
coordinatesToTSV.py $outputBase\_allCoordinates.txt $outputBase\_allCoordinates.tsv
findSuspiciousScaffolds.py data/refGenomeContigs.size $outputBase\_allCoordinates.tsv $outputBase\_highlyAlteredScaffolds.csv

#combine results into a single chromosome-level fasta file
makeFullChromFasta.sh $round $chromosomeResultsDir $outputBase

#get pbsuite stats
summarizeAssembly.py $outputBase\_chroms.fa > $outputBase\_chroms.pbSummary

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

#extract all Discovar segments, make fasta, and repeat mask
grep EI-a $outputBase\_allCoordinates.tsv|awk '$4<$5'|awk '{print $3"\t"$4"\t"$5"\tinsert\t.\t"$6}' > $outputBase\_discoInserted.bed
bedtools getfasta -s -name -fi data/DISCOVARGenome.fasta -bed $outputBase\_discoInserted.bed -fo $outputBase\_discoInserted.fa
sbatch repeatMask.slurm $outputBase\_discoInserted.fa

#get mapping quality analysis
mappingQualityAnalysis.R $overviewBedFile

#move results into the results folder.
mv $outputBase* overallResults
mv *.pdf overallResults
