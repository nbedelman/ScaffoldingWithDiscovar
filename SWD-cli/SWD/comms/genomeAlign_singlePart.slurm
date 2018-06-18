#!/bin/bash
#
# genomeAlign_last.slurm from tophat.sh example
# This script aligns the new assembly to the reference genome, and creates
# custom bed files for the next steps of the pipeline. It was created to run
# on a an assembly that has been split into several parts, and it runs the
# lastal alignment and mafToBed script on each part separately.

#SBATCH -J lastal         # A single job name for the array
#SBATCH -p serial_requeue,shared,general  # Partition
#SBATCH -n 8             # eight cores
#SBATCH -t 1-00:00         # Running time of 1 day
#SBATCH --mem 5000        # Memory request of 5 GB
#SBATCH -o outs/lastal%A_%a.out # Standard output
#SBATCH -e errs/lastal%A_%a.err # Standard error

export PATH=$PATH:/n/mallet_lab/edelman/software/kentUtils/bin
#export PATH=$PATH:/n/regal/mallet_lab/edelman/dryasAssembly/scripts

config=$1
part=$2
source $config

currentDir=$(pwd)
#create the output directories
mkdir -p $outBase
mkdir -p $outBase/sizes
mkdir -p $outBase/mafs
mkdir -p $outBase/beds
mkdir -p $outBase/net
mkdir -p $outBase/chain
mkdir -p $outBase/axt
mkdir -p $outBase/twoBit
mkdir -p $outBase/chainedAxt
mkdir -p $outBase/chainedMaf


query=$queryLoc/*.part-$part.fasta
alignOut=$outBase/mafs/$outBase.part-$part\_alns.maf
referenceSize=$outBase/sizes/$(basename $referenceAssembly).size
querySize=$outBase/sizes/$(basename $queryAssembly).size
referenceTwoBit=$outBase/twoBit/$(basename $referenceAssembly).2bit
queryTwoBit=$outBase/twoBit/$(basename $queryAssembly).2bit
axtOut=$outBase/axt/$outBase.part-$part\_alns.axt
axtSorted=$outBase/axt/$outBase.part-$part\_alns_sorted.axt
chainOut=$outBase/chain/$outBase.part-$part\_alns_sorted.chain
refNetOut=$outBase/net/$outBase.part-$part\_alns_sorted_ref.net
queryNetOut=$outBase/net/$outBase.part-$part\_alns_sorted_query.net
chainedAxtRefOut=$outBase/chainedAxt/$outBase.part-$part\_alns_sorted_chained_ref.axt
chainedAxtQuerOut=$outBase/chainedAxt/$outBase.part-$part\_alns_sorted_chained_quer.axt
chainedMafRefOut=$outBase/chainedMaf/$outBase.part-$part\_alns_sorted_chained_ref.maf
chainedMafQuerOut=$outBase/chainedMaf/$outBase.part-$part\_alns_sorted_chained_quer.maf

#find local alignments between query sequences and reference sequences
SWD.py lastAlign $referenceDB $query $lastalConfig | maf-sort.sh > $alignOut

#follow the chain and net protocol; modified from  http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto
#first, make the genome sizes file if they don't exist
if [ ! -f $referenceSize ]; then
    faSize $referenceAssembly -detailed > $referenceSize
fi
if [ ! -f $querySize ]; then
    faSize $queryAssembly -detailed > $querySize
fi

#then, convert the genomes to 2bit format if that hasn't happened yet
if [ ! -f $referenceTwoBit ]; then
    faToTwoBit $referenceAssembly $referenceTwoBit
fi
if [ ! -f $queryTwoBit ]; then
    faToTwoBit $queryAssembly $queryTwoBit
fi
#convert the last output to axt
maf-convert axt $alignOut > $axtOut
#sort the axt
axtSort $axtOut $axtSorted
#make chains
axtChain -linearGap=medium -scoreScheme=$scoreScheme $axtSorted $referenceTwoBit $queryTwoBit $chainOut
#make nets
chainNet -minScore=2000 $chainOut $referenceSize $querySize $refNetOut $queryNetOut
#convert the net back into axt
netToAxt -maxGap=300 $refNetOut $chainOut $referenceTwoBit $queryTwoBit $chainedAxtRefOut
netToAxt -maxGap=300 $queryNetOut $chainOut $referenceTwoBit $queryTwoBit $chainedAxtQuerOut
#sort and convert the chained axt back into maf
axtSort $chainedAxtRefOut  $chainedAxtRefOut.sorted.axt
axtSort $chainedAxtQuerOut  $chainedAxtQuerOut.sorted.axt
axtToMaf $chainedAxtRefOut.sorted.axt $referenceSize $querySize  $chainedMafRefOut
axtToMaf $chainedAxtQuerOut.sorted.axt $querySize $referenceSize  $chainedMafQuerOut

# ##Make a folder in which to put the bed files for each contig; then make the files
mkdir -p $outBase/beds/$part\_tmp
cd $outBase/beds/$part\_tmp

SWD.py mafToBed ../../../$chainedMafRefOut ../../../$queryAssembly

#put a copy of each bed file in the general bed folder, and concatenate into the allContigs file.
#Also put the overviews into the allOverviews file. Then get rid of the temporary folder.
for i in *.bed
do cp $i ../$i
cat $i >> $part\.bed
tail -n +2 $i >> ../$outBase\_allContigs.bed
head -n 1 $i >> ../$outBase\_allOverviews.bed
done

cd ..
rm -Rf $part\_tmp

#clean up intermediate files
cd $currentDir
rm -f $alignOut $axtOut $axtSorted $chainOut $refNetOut $queryNetOut $chainedAxtRefOut $chainedAxtQuerOut $chainedMafRefOut $chainedMafQuerOut
