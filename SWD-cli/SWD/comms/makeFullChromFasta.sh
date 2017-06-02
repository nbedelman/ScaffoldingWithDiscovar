#!/bin/bash

round=$1
inputDir=$2
outputBase=$3

for f in $inputDir/*discoOrder.fasta;
do header=$(basename $f .fasta)\_$round;
echo ">"$header >> $outputBase\_chroms.fa;
grep -v ">" $f | tr -d '\n'|awk 1 >> $outputBase\_chroms.fa;
done
formatFasta.py $outputBase\_chroms.fa
mv $outputBase\_chroms.fa\_tmp $outputBase\_chroms.fa
