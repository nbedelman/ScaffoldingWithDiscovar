#!/usr/bin/env Rscript

library(reshape2)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
genome <- args[1]

alignmentStats  <- function(f){

base=tail(strsplit(f,"/")[[1]],n=1)

#read the file
aligns=read.delim(f,sep="\t", header=F,col.names=c("Chr","ConStart","ConEnd","ID","Score","Strand","AlignStart","AlignEnd","RGB"))
aligns$Score=as.numeric(as.character(aligns$Score))


#make a data frame that has the name of each category and the number of contigs that
#fall into that category
categories=c(">5XDup","<=5XDup","2XDup","<50% Scaffold",">50% Scaffold",
             ">50% Inside",">75% Inside","All Inside","Unique Map")
summary  <- data.frame(labels=categories)
numContigs=c()
for (r in seq(200,1000,100)){
  numContigs=append(numContigs,nrow(subset(aligns,Score==r)))
}
summary=cbind(summary,numContigs)

####################### FIRST GRAPH #################################################
#Graph the overall numbers of contigs that fall into each quality category

#assign a color scheme for worst to best
colors=c("darkblue","blue", "cyan", "red","pink","orange", "yellow","green","forestgreen")
summary$labels <- factor(summary$labels,levels = categories)
# barplot(summary$numContigs,names.arg=summary$labels, col=colors,
#         main=paste(base,"Distribution of Mapping Quality"), xlab="Frequency", ylab="Quality")+
#         theme(axis.text.x = element_text(angle = 60, hjust = 1))
overviewBar <- ggplot(summary)+
  geom_bar(aes(x=labels,y=numContigs,fill=labels),col="black", stat="identity") +
  scale_fill_manual(values=colors)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(legend.title=element_blank())+
  guides(fill=FALSE)+
  xlab("Quality")+
  ylab("Number of Contigs")
overviewBar
ggsave(paste0(base,"_OverviewBarPlot.pdf"),width=11, height=8.5)

######################################################################################

####################### SECOND GRAPH #################################################
#Graph the percentage contigs that have each quality for each size bin

#define the length of an alignment as it's end position - start position
aligns$ConEnd <- as.numeric(as.character(aligns$ConEnd))
aligns$ConStart <- as.numeric(as.character(aligns$ConStart))
aligns$length=aligns$ConEnd-aligns$ConStart

#subset the alignments by length, and put the resulting data frames into a list
FiveKB=subset(aligns,length<5000)
TenKB=subset(aligns,length<10000 & length >= 5000)
FiftyKB=subset(aligns,length<50000 & length >= 10000)
HundredKB=subset(aligns,length<100000 & length >= 50000)
HundredKBPlus=subset(aligns,length>=100000)

sizes  <- list(FiveKB,TenKB,FiftyKB,HundredKB,HundredKBPlus)
sizeStrings  <- c("FiveKB","TenKB","FiftyKB","HundredKB","HundredKBPlus")

#function to get the values of each quality category for a given size category
getStats  <- function(aSubset){
  stats  <-  hist(as.numeric(as.character(aSubset$Score)), breaks=seq(100,1000,100),plot=F)$counts
  return (stats)
}

#run the above function on each size category, and put the results in a data frame.
allStats=data.frame(labels=c(">5XDup","<=5XDup","2XDup","<50% Scaffold",">50% Scaffold",">50% Inside",">75% Inside","All Inside","Unique Map"))
for (size in seq(1,length(sizes))){
  assign("x",getStats(sizes[[size]]))
  allStats=cbind(allStats,x)
}
colnames(allStats) <- c("labels",sizeStrings)

#make the absolute values into percentages
for (i in seq(2,6)){
  allStats[,i]=allStats[,i]/sum(allStats[,i])
}

#flip the data frame
n<-allStats$labels
allStats<-as.data.frame(t(allStats[,2:6]))
colnames(allStats) <- n
allStats$labels  <- factor(row.names(allStats))

#make stacked bar graph
DF2 <- melt(allStats, id.var="labels")
DF2$labels<-factor(DF2$labels,levels=c("FiveKB","TenKB","FiftyKB","HundredKB","HundredKBPlus"))
colors=c("darkblue","blue","cyan","red","pink","orange","yellow","green","forestgreen")
ggplot(DF2, aes(x = labels, y = value, fill = variable))+
  geom_bar(stat = "identity", width=.5)+ scale_fill_manual(values=colors)+
  guides(fill=FALSE)+
  xlab("Length of Contig") + ylab("Percent of Contigs") #+ ggtitle(paste(base,"Distribution of Mapping Qualities by Contig Size"))
  #theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(paste0(base,"_StackedBars.pdf"),width=11, height=8.5)

######################################################################################

####################### THIRD GRAPH ##################################################
#Make graph that shows the cumulative genomic proportion for each mapping quality, starting with best
cumLength=c()
for (i in seq(200,1000,100)){
  len=sum(subset(aligns,Score>=i)$length, na.rm=T)
  cumLength=append(cumLength,len)
}
cumLength=append(cumLength,0)
categories=append(categories,"End",after=0)
cumulScores=data.frame(x=categories,y=cumLength)

cumulScores$x<-factor(cumulScores$x,levels=c("Unique Map","All Inside",">75% Inside",">50% Inside",">50% Scaffold","<50% Scaffold","2XDup","<=5XDup",">5XDup","End"))

colors2a=append(rev(colors),"darkblue")
ggplot(data=cumulScores, aes(x=x, y=y, group=1)) +
  geom_line(aes(color=x), size=3)+ scale_color_manual(values=colors2a) +
  xlab("Map Quality") + ylab("Cumulative Length") +
  guides(color=FALSE)+
  #geom_hline(aes(yintercept=2.76e+08),linetype="dashed") +
  xlab("Quality")+
  ylab("Length")
  #labs(list(title=paste0(base," cumulative length by quality"),x = "quality", y="length"))#+
  #theme(legend.position="none")
ggsave(paste0(base,"_LineGraph.pdf"),width = 11,height = 8.5)
######################################################################################
}

#args=(commandArgs(TRUE))
#for (i in 1:length(args)){
#  eval(parse(text=args[[i]]))
#}

alignmentStats(genome)
