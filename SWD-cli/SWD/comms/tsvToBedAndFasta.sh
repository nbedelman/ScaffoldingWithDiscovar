#This script takes the new order tsv file and creates a bed file that can be used by
#bedtools getFasta to create the newly ordered fasta file for the whole genome.
#Then, it makes the fasta file.

tsvFile=$1
outputBase=$2
genome=$3
bedFile=$outputBase\.bed
fastaFile=$outputBase\.fa

#First, add 'Chr' to the chromosome number line
#this only makes sense for Hmel
#nl -s Chr $tsvFile| cut -c7- |\

#Then, get rid of the header
tail -n +2 $tsvFile|\

#then, only output the relevant lines
awk -v OFS='\t' '{print $3, $4, $5, $1,"0",$6}' > $bedFile

#next, use bedtools to convert the bed to a fasta. -s uses the strand info.
bedtools getfasta -s -name -fi $genome -bed $bedFile -fo $fastaFile\_tmp

#finally, get rid of the extraneous fasta headers and re-format
python removeExcessHeaders.py $fastaFile\_tmp > $fastaFile
python formatFasta.py $fastaFile
mv $fastaFile\_tmp $fastaFile
