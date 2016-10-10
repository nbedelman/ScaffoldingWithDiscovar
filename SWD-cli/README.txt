# ScaffoldingWithDiscovar

Set of scripts to use DISCOVAR data to order and orient scaffolds in a genome. The genome must have some scaffolds assigned to chromosomes.


Required external software:
LAST (http://last.cbrc.jp/)
Bedtools (http://bedtools.readthedocs.org/en/latest/)
bioPython (https://github.com/biopython/biopython)

This program is focused on the "low-hanging fruit" of genomic scaffolding. It takes a de-novo assembled genome and a reference genome, and simply sees
if there are any instances where the new genome assembly has a scaffold that would connect two scaffolds of the reference. It only considers reference
scaffolds "bridgeable" if they have already been assigned to the same chromosome, though it would not be terribly difficult to adjust that. It does make 
for a stronger filter, though, and we can be more confident in the result.


##PREPARING THE DATA###

Start with  the output from a novel genome assembly (fasta), a fasta file of reference genome scaffolds (fasta), and a map file for the reference genome (bed). 

Reference Genome Map:
This should be in a bed format (see https://genome.ucsc.edu/FAQ/FAQformat#format1)

example:
<chrom> <start> <stop>  <name>  <score> <strand>
Chr1	1	73148	Hmel201001	500	+
Chr1	73149	73248	Ns	500	+
Chr1	73249	2694103	Hmel201002	1000	+

*In this example, when the chromosome is built, each scaffold is separated by 100 Ns. This is accomplished by including an entry called Ns and including it between each scaffold. This is not necessary, and separation of scaffolds by Ns can be done in other ways as well. This particular method is not necessary for SWD to work.

Reference Genome Sequence:
This should be a fasta file with each scaffold as an entry. They should be in "positive" orientation.

Novel Genome Assembly:
Only use contigs that have a length >= 1000 bp (removeSmallScaffolds.py)

###ALIGNING SEQUENCES###

The tool SWD lastAlign uses the program LAST, and SWD was only tested with LAST. However, theoretically any aligner should be fine,
as long as you get the output in MAF format. 

***The novel assembly must be mapped to full chromosomes. SWD genomeBuild combines scaffolds into chromosomes using the provided map file before creating a genome database. If you are using another alignment program, please keep this in mind and make sure to create a chromosomal fasta file prior to alignment.***

In this section, I will assume you are using LAST. If that's not the case, use another 
aligner and skip to the next section.

To use LAST, you will have to create a reference genome database with lastdb and then align a query to the reference with lastal. 
To hopefully make this easier, I have created SWD commands for these processes (SWD genomeBuild and SWD lastAlign) with config files that you can edit. They simply run LAST, but nicely within the SWD suite. For LAST documentation, see http://last.cbrc.jp/

If the genome assembly is large (likely), first split it for better parallelization (fasta-splitter.pl). Then, align each part to genome with SWD lastAlign.

The output of this step will be maf files for each of the query sequences you aligned.

###CONVERTING MAF TO BED AND FINDING RELEVANT SCAFFOLDS###

Once you have the MAF alignments, convert them into BED files using SWD mafToBed. Each new assembly scaffold will have a bed file containing the top
hits for that sequence, as well as a custom score.
Then, pull out the new assembly scaffolds that might be informative for re-scaffolding the reference with SWD overviewToIntersect
***at this point, can check quality of alignment with mappingQualityAnalysis.R

The output of this step will be a directory called overlapBed which contains the bed files of interest.

###RE-ORDER THE REFERENCE GENOME###

