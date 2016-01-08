# ScaffoldingWithDiscovar

Set of scripts to use DISCOVAR data to order and orient scaffolds in a genome. The genome must have some scaffolds assigned to chromosomes.


Required external software:
LAST (http://last.cbrc.jp/)
Bedtools (http://bedtools.readthedocs.org/en/latest/)
bioPython (https://github.com/biopython/biopython)

Start with  the output from DISCOVAR de novo (fasta), a scaffolded genome (fasta), and a map file for that genome (agp)

1) only take those DISCOVAR contigs that have a length >= 1000 bp (removeSmallScaffolds.py)
2) use LAST to align the DISCOVAR data to the scaffolded genome
  a) if the DISCOVAR data is large (likely), first split it for better parallelization (fasta-splitter.pl)
  b) index the genome with lastdb
  c) align each part to genome with lastal
  d) convert the resulting alignments in maf format to bed format. Only take the best alignments
      - during this step, make both individual bed files for each contig AND the custom "overview" bed file
  ***b-d all done with genomeAlign_last_array.slurm
  ***at this point, can check quality of alignment with mappingQualityAnalysis.R
3) find contigs that overlap more than one reference scaffold (overviewToIntersects.sh)
4) use the overlapping alignments to re-order the reference scaffolds (reOrderChroms.py)
