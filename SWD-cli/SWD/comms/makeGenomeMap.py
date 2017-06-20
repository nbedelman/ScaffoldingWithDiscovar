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



def makeInexactGenomeMap(genome,Nlength,namePrefix="name"):
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
            if len(seqChunks[chunk]) > 1:
                mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(seqChunks[chunk]),chromName+"_"+str(contig),"500","+"))
            pos=pos+len(seqChunks[chunk])
            if len(Nchunks[chunk]) > 1:
                mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(Nchunks[chunk]),"Ns_"+str(nNum),"500","+"))
            pos=pos+len(Nchunks[chunk])
            contig+=1
            nNum+=1
        if len(seqChunks[-1]) > 1:
            mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(seqChunks[-1]),chromName+"_"+str(contig),"500","+"))
    mapFile.close()

def makeExactGenomeMap(genome,Nlength,namePrefix="name"):
    chromNum=0
    mapFile=open(genome+".bed","w")
    g=SeqIO.parse(genome,"fasta")
    Npattern="[ACTG]"+"N"*Nlength+"[ACTG]"
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
        #write the first entry. If it is the ONLY entry, write the whole thing without adding one to the end.
        if len(seqChunks)==1:
            mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(seqChunks[0]),chromName+"_"+str(contig),"500","+"))
            pass
        else:
            mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(seqChunks[0])+1,chromName+"_"+str(contig),"500","+"))
            pos=pos+len(seqChunks[0])+1
            contig+=1

        for chunk in range(len(Nchunks)-1):
            if len(Nchunks[chunk]) > 1:
                mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(Nchunks[chunk])-2,"Ns_"+str(nNum),"500","+"))
                pos=pos+len(Nchunks[chunk])-2
            if len(seqChunks[chunk+1]) > 1:
                mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(seqChunks[chunk+1])+2,chromName+"_"+str(contig),"500","+"))
                pos=pos+len(seqChunks[chunk+1])+2
            contig+=1
            nNum+=1

        #write the last entries
        try:
            mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(Nchunks[-1])-2,"Ns_"+str(nNum),"500","+"))
            pos=pos+len(Nchunks[-1])-2
            mapFile.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (chromName, str(pos), pos+len(seqChunks[-1])+1,chromName+"_"+str(contig),"500","+"))
            nNum+=1
        except IndexError:
            #If there are no stretches of Ns that have the exact length desired, the whole chromosome will be written as one entry above.
            pass
    mapFile.close()

genome=sys.argv[1]
Nlength=int(sys.argv[2])
namePrefix=sys.argv[3]
try:
    exact=sys.argv[4]
except IndexError:
    exact=False

if exact:
    makeExactGenomeMap(genome,Nlength,namePrefix)
else:
    makeInexactGenomeMap(genome,Nlength,namePrefix)
