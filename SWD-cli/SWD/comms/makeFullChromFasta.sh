#!/bin/bash

round=$1
outputFile=$2

for f in *discoOrder.fasta
do header=$(basename $f .fasta)\_$round
echo $header >> $outputFile
grep -v ">" $f >> $outputFile
done
