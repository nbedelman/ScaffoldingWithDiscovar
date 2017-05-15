#take an overview bed file, and filter it, keeping only good alignments that overlap >1 reference scaffold.
#In this case, "good" means non-repetitive sequence and overlapping multiple reference scaffolds from the same
#chromosome.

overviewBedFile=$1
agpBedFile=$2
bedDirectory=$3

mkdir -p overlapBeds

awk '$2>0' $overviewBedFile | awk '$5~/600|700|800/' | sed 's/ ./\t/g' > $overviewBedFile\_goodScores.bed
sed 's/ ./\t/g' $agpBedFile > $agpBedFile\_formatted.bed
bedtools intersect -c -a $overviewBedFile\_goodScores.bed -b $agpBedFile\_formatted.bed | awk '$10 > 1' | cut -f 4 > overlaps.txt

while read line;
do bedFile=$line.bed
cp $bedDirectory/$bedFile overlapBeds/
done < overlaps.txt
