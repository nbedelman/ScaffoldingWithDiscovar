#Take an AGP file, output a bed file that can be used with bedTools to give a fasta file.
from Bio.Seq import Seq
from Bio import SeqIO

def ChrToHmel(chromo):
    '''Takes a chromosome ID in format ChrX or ChrXX
    returns chromosomeID in format Hmel2XX'''
    num=str(chromo[3:]).strip("_unmapped")
    if len(num) ==1:
        return "Hmel20"+num
    elif len(num)==2:
        return "Hmel2"+num

class agpEntry(object):
    def __init__(self,entry,color=False):
        colorList=["150,150,150","204,204,204","225,225,225"]
        atts=entry.split("\t")
        self.chrom=ChrToHmel(atts[0])
        self.start=atts[1]
        self.end=atts[2]
        self.num=atts[3]
        self.DN=atts[4]
        self.scaf=atts[5]
        self.oneFrag=atts[6]
        if self.DN == "D":
            self.length=atts[7]
            self.strand=atts[8].strip("\n")
        elif self.DN =="N":
            self.length=100
            self.strand="+"
        try:
            self.status=atts[9]
            self.phys=atts[10]
        except IndexError:
            self.status="NA"
            self.phys="NA"
        if type(color)==int:
            self.color=colorList[color]
    def getChrom(self):
        return self.chrom
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getNum(self):
        return self.num
    def getDN(self):
        return self.DN
    def getScaf(self):
        return self.scaf
    def getOneFrag(self):
        return self.oneFrag
    def getLength(self):
        return self.length
    def getStrand(self):
        return self.strand
    def getStatus(self):
        return self.status
    def getPhys(self):
        return self.phys
    def getColor(self):
        return self.color
    def writeBed(self,scaffs=False, color=False):
        score=0
        if "anchored" in self.getStatus():
            score=1000
        elif "orient" in self.getStatus():
            score=500
        outList=[]
        if not scaffs:
            outList+=[self.getChrom(),str(self.getStart()),str(self.getEnd()),self.getScaf(),score,self.getStrand()]
        else:
            outList+=[self.getChrom(),"1",str(self.getLength()),self.getScaf(),score,self.getStrand()]
        if color:
            outList+=[outList[1],outList[2],self.getColor()]
        toWrite=''
        for item in outList:
            toWrite += str(item)+"\t"
        toWrite+="\n"
        if len(outList)==6:
            toWrite+=self.getChrom()+"\t"+"1"+"\t"+"100"+"\t"+"Ns"+"\t"+"500"+"\t"+"+"+"\n"
        elif len(outList)==9:
            toWrite+=self.getChrom()+"\t"+"1"+"\t"+"100"+"\t"+"Ns"+"\t"+"500"+"\t"+"+"+"\n"+"1"+"\t"+"100"+"\t"+"100,100,100"+"\n"
        return toWrite

def readAGP(agpFile, color=False):
    '''take an agp file
    return a list of agpEntries, in the order they are in the file'''
    agpList=[]
    with open(agpFile,"r") as f:
        count=0
        for line in f:
            if color:
                agpList.append(agpEntry(line, color=count%3))
                count+=1
            else:
                agpList.append(agpEntry(line))
    return agpList

def addUnmappedChroms(bedFile,fastaFile, color=False, scaffs=False):
    with open(bedFile, "r") as b:
            allScaffs=[line.split()[3] for line in b]
    b.close()
    f = SeqIO.parse(open(fastaFile, "r"), "fasta")
    start=1
    end=1
    count=0
    colorList=["150,150,150","204,204,204","225,225,225"]
    with open(bedFile, "a+") as o:
        for record in f:
            if record.id not in allScaffs:
                scaf=record.id
                length=len(record)
                strand="+"
                toWrite=[]
                if scaffs:
                    toWrite+=["Hmel200","1",str(length),scaf,"0",strand]
                else:
                    end=start+length
                    toWrite+=["Hmel200",start,end,scaf,"0",strand]
                    start=end+1
                if color:
                    toWrite+=[start,end,colorList[count%3]]
                    count+=1
                output=''
                for i in toWrite:
                    output+=str(i)+"\t"
                output += "\n"
                o.write(output)
    o.close()
    f.close()

def writeBedForFasta(agpFile,fastaFile, outFile):
    scafList=readAGP(agpFile, color=True)
    with open(outFile,"w") as f:
        for entry in scafList:
            if entry.getDN() == "D":
                f.write(entry.writeBed(color=False,scaffs=True))
    f.close()
    addUnmappedChroms(outFile, fastaFile, color=False, scaffs=True)

def convertScaffsToChroms(bedFile):
    names=bedFile[:-4]
    f=open(bedFile, "r")
    start=1
    end=1
    currentChrom=''
    with open(names+"_chroms.bed", "w") as o:
        for line in f:
            atts=line.split("\t")
            chrom=atts[0]
            if chrom==currentChrom:
                length=int(atts[2])
                end=start+(length-1)
                outString=str(atts[0])+"\t"+str(start)+"\t" +str(end)+"\t"+str(atts[3]) + "\t"+str(atts[4]) + "\t"+str(atts[5].strip("\n"))+"\n"
                o.write(outString)
                start=end+1
            else:
                currentChrom=chrom
                start=1
                length=int(atts[2])
                end=start+(length-1)
                outString=str(atts[0])+"\t"+str(start)+"\t" +str(end)+"\t"+str(atts[3]) + "\t"+str(atts[4]) + "\t"+str(atts[5].strip("\n")) + "\n"
                o.write(outString)
                start=end+1
        o.close
    f.close()
