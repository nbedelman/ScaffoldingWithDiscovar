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

Start with  the output from a novel genome assembly (fasta), a scaffolded reference genome (fasta), and a map file for the reference genome (agp)

Reference Genome Map:
This should be an AGP file in the format:
<chromosomeName>	<start>	<end>	<number>	<D/N>	<scaffoldName>	<1>	<length>	<strand>	<optional additional columns>


Reference Genome Sequence:
This program uses as input a reference genome in fasta format, where each sequence is an entire chromosome. If your genome is 
a fasta file that has a different entry for each scaffold, you will need to combine the scaffolds into chromosome. To do this, 
you can use SWD agpToFasta (uses the scripts agpToBedForFasta.py, and then bedtools getfasta).

Novel Genome Assembly:
Only use contigs that have a length >= 1000 bp (removeSmallScaffolds.py)

###ALIGNING SEQUENCES###

The tool SWD lastAlign uses the program LAST, and SWD was only tested with LAST. However, theoretically any aligner should be fine,
as long as you get the output in MAF format. In this section, I will assume you are using LAST. If that's not the case, use another 
aligner and skip to the next section.

To use LAST, you will have to create a reference genome database with lastdb and then align a query to the reference with lastal. 
To hopefully make this easier, I have created SWD commands for these processes (SWD lastDB and SWD lastAlign) with config files that
you can edit. They simply run LAST, but nicely within the SWD suite. For LAST documentation, see http://last.cbrc.jp/

If the genome assembly is large (likely), first split it for better parallelization (fasta-splitter.pl). Then, align each part to genome with lastal.

The output of this step will be maf files for each of the query sequences you aligned.

###CONVERTING MAF TO BED AND FINDING RELEVANT SCAFFOLDS###

Once you have the MAF alignments, convert them into BED files using SWD mafToBed. Each new assembly scaffold will have a bed file containing the top
hits for that sequence, as well as a custom score.
Then, pull out the new assembly scaffolds that might be informative for re-scaffolding the reference with SWD overviewToIntersect
***at this point, can check quality of alignment with mappingQualityAnalysis.R

The output of this step will be a directory called overlapBed which contains the bed files of interest.

###RE-ORDER THE REFERENCE GENOME###

