#take an overview bed file, and filter it, keeping only good alignments that overlap >1 

discoBedFile=$1
agpBedFile=$2

mkdir overlapBeds

grep -v -E "\-\d" $discoBedFile | grep -e "\t600\t" -e "\t700\t" -e "\t800\t"\
 | bedtools intersect -c -a - -b $agpBedFile \
 | awk '$10 > 1' | cut -f 4 > overlaps.txt

while read line;
do bedFile=$line.bed
cp fullOverlaps/$bedFile overlapBeds/
done < overlaps.txt

