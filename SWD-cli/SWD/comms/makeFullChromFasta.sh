#!/bin/env python

import sys
import os
from Bio import SeqIO

round=sys.argv[1]
inputDir=sys.argv[2]
outputBase=sys.argv[3]

fullChroms=[]
for subdir, dirs, files in os.walk(inputDir):
    for file in files:
        if "discoOrder.fasta" in file:
            header=file.split(".")[0]+"_"+str(round)
            prevNs=False
            chrom=SeqIO.parse(inputDir+"/"+file,"fasta")
            thisChrom=SeqIO.SeqRecord(seq="", id=header)
            for record in chrom:
                if "Ns_" in record.id:
                    prevNs=True
                    thisChrom.seq+=record.seq
                else:
                    if not prevNs:
                        thisChrom.seq+="N"*500
                    thisChrom.seq+=record.seq
                    prevNs=False
            fullChroms.append(thisChrom)
SeqIO.write(fullChroms,outputBase+"_chroms.fa", "fasta")
