#!/usr/bin/env python

from Bio import SeqIO
import re
import sys

# ordered="D.iulia_2.0.fasta"
# unordered="D.iulia_unordered.fasta"
# ungrouped="D.iulia_unclustered.fasta"

def combineFastas(ordered,unordered,ungrouped):
    o=SeqIO.parse(ordered,"fasta")
    u=SeqIO.index(unordered,"fasta")
    g=SeqIO.parse(ungrouped,"fasta")
    out=open("combinedGenome.fasta","w")
    for record in o:
        newRecord=SeqIO.SeqRecord(seq=record.seq,id=record.id.split("__")[0],name=record.id.split("__")[0])
        grp="group"+record.id.split("group")[1].split("_")[0]
        for k in u.keys():
            if k.endswith(grp):
                newRecord.seq+="N"*50000
                newRecord.seq+=u[k].seq
        SeqIO.write([newRecord,],out,"fasta")
    unClust=SeqIO.SeqRecord(seq="",id="Lachesis_group31",name="unclustered")
    for r in g:
        unClust.seq+=r.seq
        unClust.seq+="N"*50000
    SeqIO.write([unClust,],out,"fasta")
    o.close()
    u.close()
    g.close()
    out.close()

#combineFastas(ordered,unordered,ungrouped)



genome=sys.argv[1]
Nlength=int(sys.argv[2])
namePrefix=sys.argv[3]

def makeGenomeMap(genome,Nlength,namePrefix="name"):
    chromNum=0
    mapFile=open(genome+".bed","w")
    g=SeqIO.parse(genome,"fasta")
    Npattern="N"*Nlength+"*"
    nNum=0
    for chrom in g:
        pos=0
        contig=0
        if namePrefix=="name":
            chromName=chrom.id
        else:
            chromName=namePrefix+str(chromNum)
        chromNum+=1
        stringChrom=str(chrom.seq)
        Nchunks=re.findall(Npattern,stringChrom)
        seqChunks=re.split(Npattern,stringChrom)
        for chunk in range(len(Nchunks)):
            if pos+len(seqChunks[chunk]) > pos+1:
                mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos+1), pos+len(seqChunks[chunk]),chromName+"_"+str(contig),"500","+"))
            pos=pos+len(seqChunks[chunk])
            if pos+len(Nchunks[chunk]) > pos+1:
                mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos+1), pos+len(Nchunks[chunk]),"Ns_"+str(nNum),"500","+"))
            pos=pos+len(Nchunks[chunk])
            contig+=1
            nNum+=1
        if pos+len(seqChunks[-1]) > pos+1:
            mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos+1), pos+len(seqChunks[-1]),chromName+"_"+str(contig),"500","+"))
    mapFile.close()

makeGenomeMap(genome,Nlength,namePrefix)
