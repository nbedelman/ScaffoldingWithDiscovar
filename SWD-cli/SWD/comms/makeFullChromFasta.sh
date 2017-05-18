#!/bin/bash

round=$1
outputFile=$2

for f in *discoOrder.fasta
do header=$(basename $f .fasta)\_$round
echo ">"$header >> $outputFile
grep -v ">" $f | tr -d '\n'|awk 1 >> $outputFile
done
python ~/Documents/Mallet_Lab/referenceScaffolding/ScaffoldingWithDiscovar/SWD-cli/SWD/comms/formatFasta.py $outputFile
#mv $outputFile\_tmp $outputFile