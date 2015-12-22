#take an overview bed file, and filter it by getting rid of negative values

grep -v -E "\-\d" h_melpomene_allOverviews.bed | grep -e "\t600\t" -e "\t700\t" -e "\t800\t"\
 | bedtools intersect -c -a - -b agpToBed_chroms.bed \
 > multiIntersects.bed
