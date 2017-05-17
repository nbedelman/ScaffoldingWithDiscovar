#formatFasta.py

#this seems silly, but easiest way I know to make sure the fasta is properly formatted

from Bio import SeqIO
import sys

inFile=sys.argv[1]

genome=SeqIO.parse(inFile, "fasta")
outputFile=inFile+"_tmp"
chroms=[]
for record in genome:
    chroms.append(record)
SeqIO.write(chroms,outputFile, "fasta")