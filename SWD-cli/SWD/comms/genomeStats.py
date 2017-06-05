#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
import csv
import os
import sys

def getGenomeStats(genomeFile, format):
    lengths=[]
    cumuLengths=[]
    ge1KB=[]
    ge5KB=[]
    if format=="fasta":
        genome=SeqIO.parse(open(genomeFile, "r"), "fasta")
        for record in genome:
            lengths.append(len(record))
            if len(record)>=1000:
                ge1KB.append(len(record))
                if len(record)>=5000:
                    ge5KB.append(len(record))
    elif format=="bed":
        genome=open(genomeFile,"r")
        for entry in genome:
            atts=entry.split("\t")
            length=int(atts[2])-int(atts[1])+1
            lengths.append(length)
            if length >= 1000:
                ge1KB.append(length)
                if length>=5000:
                    ge5KB.append(length)
    total=0
    lengths=sorted(lengths)
    for i in lengths:
        cumuLengths.append(i+total)
        total+=i
    totSize=cumuLengths[-1]
    n50Index=next(i for i, x in enumerate(cumuLengths) if x > (totSize*0.5))
    n90Index=next(i for i, x in enumerate(cumuLengths) if x > (totSize*0.9))
    n95Index=next(i for i, x in enumerate(cumuLengths) if x > (totSize*0.95))
    n50=lengths[n50Index]
    n90=lengths[n90Index]
    n95=lengths[n95Index]
    medianSize=getMedian(lengths)
    numSeqs=len(lengths)
    minSize=lengths[0]
    meanSize=float(sum(lengths))/len(lengths)
    maxSize=lengths[-1]
    firstQ=getMedian(lengths[:int(len(lengths)/2)])
    thirdQ=getMedian(lengths[int(len(lengths)/2):])
    ge1KBseqs=len(ge1KB)
    ge1KBpct=float(sum(ge1KB))/sum(lengths)
    ge5KBseqs=len(ge5KB)
    ge5KBpct=float(sum(ge5KB))/sum(lengths)
    print '''%s,%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t''' % (genomeFile.split("/")[-1],numSeqs,\
    minSize,firstQ,medianSize,meanSize,thirdQ,maxSize,totSize,n50,n90,n95,ge1KBseqs,ge1KBpct,ge5KBseqs,ge5KBpct)
    return [genomeFile.split("/")[-1],numSeqs,\
    minSize,firstQ,medianSize,meanSize,thirdQ,maxSize,totSize,n50,n90,n95,ge1KBseqs,ge1KBpct,ge5KBseqs,ge5KBpct]

def getMedian(aList):
    if len(aList)%2==1:
        medianSpot=(len(aList)+1)/2
        median=aList[medianSpot]
    elif len(aList)%2 == 0:
        midOne=len(aList)/2
        midTwo=midOne+1
        median=(aList[midOne]+aList[midTwo])/2
    return median

# allGenomeStats=[]
# rootdir="/Users/nbedelman/Documents/Mallet_Lab/18Genomes/Genomic_Analysis/clippedGenomes"
# for subdir, dirs, files in os.walk(rootdir):
#     for file in files:
#         if '.fasta' in file:
#             allGenomeStats.append(getGenomeStats(os.path.join(subdir,file)))



def write_csv(file, asequence, header=None):
    fp =  open(file, 'w')
    a = csv.writer(fp, delimiter=',')
    if header:
        a.writerows(header)
    a.writerows(asequence)
    fp.close()

stats=getGenomeStats(sys.argv[1], sys.argv[2])
write_csv(sys.argv[3],[stats,], header=[["genomeName","numSeqs",\
"minSize","firstQ","medianSize","meanSize","thirdQ","maxSize","totSize","n50","n90","n95","ge1KBseqs","ge1KBpct","ge5KBseqs","ge5KBpct"]])
