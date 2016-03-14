library(ggplot2)

envelopers=read.csv("envelopedLengths.csv", header=FALSE, col.names = "enveloped_sizes")
hist(envelopers$enveloped_sizes, breaks = 100, main="size frequency of enveloped scaffolds", xlab = "size")
sum(envelopers$enveloped_sizes)
median(envelopers$enveloped_sizes)
nrow(envelopers)

ggplot(data=envelopers, aes(enveloped_sizes))+
  geom_histogram(binwidth=.1)+
  scale_x_continuous(trans="log", limits=c(1000,1000000), breaks=c(1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000))+
  labs(title="size frequency of enveloped scaffolds", x="size (log scale)")

original=read.csv("genomeScafLengths.csv", header=FALSE, col.names = "scaf_sizes")
ggplot()+
  geom_histogram(data=original, aes(scaf_sizes),binwidth=.1)+
  geom_histogram(data=unique, aes(V8),binwidth=.1, color="blue")+
  scale_x_continuous(trans="log", limits=c(200,10000000), 
                     breaks=c(200,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000, 5000000, 10000000 ))+
  labs(title="size frequency of scaffolds", x="size (log scale)")


unique=read.csv("uniqueOverlappers.csv", header=FALSE)
