#take an overview bed file, and filter it, keeping only good alignments that overlap >1 reference scaffold.
#In this case, "good" means non-repetitive sequence and overlapping multiple reference scaffolds from the same
#chromosome.

overviewBedFile=$1
agpBedFile=$2
bedDirectory=$3

mkdir overlapBeds

grep -v -E "\-\d" $overviewBedFile | grep -e "\t600\t" -e "\t700\t" -e "\t800\t"\
 | bedtools intersect -c -a - -b $agpBedFile \
 | awk '$10 > 1' | cut -f 4 > overlaps.txt

while read line;
do bedFile=$line.bed
cp $bedDirectory/$bedFile overlapBeds/
done < overlaps.txt

