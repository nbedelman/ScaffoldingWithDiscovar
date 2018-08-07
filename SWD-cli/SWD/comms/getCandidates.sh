#take an overview bed file, and filter it, keeping only good alignments that overlap >1 reference scaffold.
#In this case, "good" means non-repetitive sequence and overlapping multiple reference scaffolds from the same
#chromosome.

overviewBedFile=$1
bedDirectory=$2

mkdir -p candidateBeds

awk '$2>0' $overviewBedFile | awk '$5~/600|700|800/' | sed 's/ ./\t/g' > $overviewBedFile\_goodScores.bed
cut -f 4 $overviewBedFile\_goodScores.bed > candidates.txt

while read line;
do bedFile=$line.bed
directions=$(awk '{print $6}' $bedDirectory/$bedFile|sort|uniq|wc -l|awk '{print $1}')
echo $directions
if [ $directions == 2 ];then
cp $bedDirectory/$bedFile candidateBeds/
fi
done < candidates.txt
