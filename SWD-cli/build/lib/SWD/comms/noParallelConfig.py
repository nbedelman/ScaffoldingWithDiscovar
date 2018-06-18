#!/usr/bin/env python


#This file allows the user to specify options for lastal. Each option is explained via its documentation entry from http://last.cbrc.jp/doc/lastal.html and then assigned a variable name. Simply change the value of the variable to specify an option.

####################
###OPTIONS I USED###
####################

#Maximum multiplicity for initial matches. Each initial match is lengthened until it occurs at most this many times in the reference.
#If the reference was split into volumes by lastdb, then lastal uses one volume at a time. The maximum multiplicity then applies to each volume, not the whole reference. This is why voluming changes the results.
#-m multiplicity
multiplicity=100

#Minimum alignment score. (If you do gapless alignment with option -j1, then -d and -e mean the same thing. If you set both, -e will prevail.)
#Initial-match options
#-e Score
minAlignScore=500

#Specify lowercase-marking of repeats, by two digits (e.g. "01"), with the following meanings.
#First digit:
#0. Convert the input sequences to uppercase while reading them.
#1. Keep any lowercase in the input sequences.
#Second digit:
#0. Do not check for simple repeats.
#1. Convert simple repeats (e.g. cacacacacacacacac) to lowercase.
#2. Convert simple repeats, within AT-rich DNA, to lowercase.
#-R Digits
repeatHandling="01"

#Specify a match/mismatch score matrix. Options -r and -q will be ignored. The built-in matrices are described in last-matrices.html.
#Any other NAME is assumed to be a file name. For an example of the format, see the matrix files in the data directory. Any letters that aren't in the matrix will get the lowest score in the matrix when aligned to anything. Asymmetric scores are allowed: query letters correspond to columns and reference letters correspond to rows. Other options can be specified on lines starting with "#last", but command line options override them.
#-p Name
scoreMatrix="ATMAP"

#####################
###E-VALUE OPTIONS###
#####################

#Report alignments that are expected by chance at most once per LENGTH query letters. This option only affects the default value of -E, so if you specify -E then -D has no effect.
#-D Length
alignExpect=None

#Maximum EG2 (expected alignments per square giga). This option only affects the default value of -e, so if you specify -e then -E has no effect.
#-E Threshold
maxEG2=None

###################
###SCORE OPTIONS###
###################

#Match score
#-r Score
matchScore=None

#Mismatch cost
#-q Cost
misMatch=None

#Gap Existence Cost
#-a cost
gapExistCost=None

#Gap Extension Cost
#-b cost
gapExtendCost=None

#Insertion existence cost. This refers to insertions in the query relative to the reference. If -A is not used, the insertion existence cost will equal the deletion existence cost, which is set by -a.
#-A Cost
insertExistCost=None

#Insertion extension cost.
#-B Cost
insertExtendCost=None

#-c Cost
afflineGap=None

#Maximum score drop for gapped alignments. Gapped alignments are forbidden from having any internal region with score < -DROP. This serves two purposes: accuracy (avoid spurious internal regions in alignments) and speed (the smaller the faster).
#-x Drop
maxScoreDropGap=None

#Maximum score drop for gapless alignments.
#-y Drop
maxScoreDropGapless=None

#Maximum score drop for final gapped alignments. Setting z different from x causes lastal to extend gapless alignments twice: first with a maximum score drop of x, and then with a (presumably higher) maximum score drop of z.
#-z Drop
maxScoreDropGapFinal=None

#Minimum score for gapless alignments.
#-d Drop
minScoreGap=None

#Minimum length for initial matches. Length means the number of letters spanned by the match.
#-l length
minInitMatchLength=None

#Maximum length for initial matches.
#-L length
maxInitMatchLength=None

#Look for initial matches starting only at every STEP-th position in each query (positions 0, STEP, 2xSTEP, etc). This makes lastal faster but less sensitive.
#-k Step
searchStep=None

#Look for initial matches starting only at query positions that are "minimum" in any window of SIZE consecutive positions (see lastdb.html). By default, this parameter takes the same value as was used for lastdb -W.
#-W Size
minQuerySize=None

###########################
###MISCELLANEOUS OPTIONS###
###########################

#Specify which query strand should be used: 0 means reverse only, 1 means forward only, and 2 means both.
#-s Strand
queryStrand=None

#Specify which DNA strand the score matrix applies to. This matters only for unusual matrices that lack strand symmetry (e.g. if the a:g score differs from the t:c score). 0 means that the matrix applies to the forward strand of the reference aligned to either strand of the query. 1 means that the matrix applies to either strand of the reference aligned to the forward strand of the query.
#-S Number
dnaStrand=None

#Omit any alignment whose query range lies in LIMIT or more other alignments with higher score (and on the same strand). This is a useful way to get just the top few hits to each part of each query (P Berman et al. 2000, J Comput Biol 7:293-302).
#-K Limit
numOverlapQuery=None

#Before extending gapped alignments, discard any gapless alignment whose query range lies in LIMIT or more others (for the same strand and volume) with higher score-per-length. This can reduce run time and output size (MC Frith & R Kawaguchi 2015, Genome Biol 16:106).
#-C Limit
numOverlapExtend=None

#Divide the work between this number of threads running in parallel. 0 means use as many threads as your computer claims it can handle simultaneously. Single query sequences are not divided between threads, so you need multiple queries per batch for this option to take effect.
#-P Threads
numThreads=None

#Search queries in batches of at most this many bytes. If a single sequence exceeds this amount, however, it is not split. You can use suffixes K, M, and G to specify KibiBytes, MebiBytes, and GibiBytes. This option has no effect on the results (apart from their order).
#If the reference was split into volumes by lastdb, then each volume will be read into memory once per query batch.
#-i Bytes
batchSize=None

#Find minimum-difference alignments, which is faster but cruder. This treats all matches the same, and minimizes the number of differences (mismatches plus gaps).
#Any substitution score matrix will be ignored. The substitution scores are set by the match score (r) and the mismatch cost (q).
#The gap cost parameters will be ignored. The gap existence cost will be 0 and the gap extension cost will be q + r/2.
#The match score (r) must be an even number.
#Any sequence quality data (e.g. fastq) will be ignored.
#-M
minDiffAligns=None

#Type of alignment: 0 means "local alignment" and 1 means "overlap alignment". Local alignments can end anywhere in the middle or at the ends of the sequences. Overlap alignments must extend to the left until they hit the end of a sequence (either query or reference), and to the right until they hit the end of a sequence.
#Warning: it's often a bad idea to use -T1. This setting does not change the maximum score drop allowed inside alignments, so if an alignment cannot be extended to the end of a sequence without exceeding this drop, it will be discarded.
#-T Number
alignType=None

#Maximum number of gapless alignments per query position. When lastal extends gapless alignments from initial matches that start at one query position, if it gets COUNT successful extensions, it skips any remaining initial matches starting at that position.
#-n Count
alignPerQueryPosition=None

#Specify treatment of lowercase letters when extending alignments:
#0. Mask them for neither gapless nor gapped extensions.
#1. Mask them for gapless but not gapped extensions.
#2. Mask them for gapless but not gapped extensions, and then discard alignments that lack any segment with score >= e when lowercase is masked.
#3. Mask them for gapless and gapped extensions.
#"Mask" means change their match/mismatch scores to min(unmasked score, 0). This option does not affect treatment of lowercase for initial matches.
#-u NUMBER
treatLower=None

#This option is a kludge to avoid catastrophic time and memory usage when self-comparing a large sequence. If the sequence contains a tandem repeat, we may get a gapless alignment that is slightly offset from the main self-alignment. In that case, the gapped extension might "discover" the main self-alignment and extend over the entire length of the sequence.
#To avoid this problem, gapped alignments are not triggered from any gapless alignment that:
#is contained, in both sequences, in the "core" of another alignment
#has start coordinates offset by DISTANCE or less relative to this core
#Use -w0 to turn this off.
#-w DISTANCE
timeKludge=None

#	Use an alternative genetic code in the specified file. For an example of the format, see vertebrateMito.gc in the examples directory. By default, the standard genetic code is used. This option has no effect unless DNA-versus-protein alignment is selected with option -F.
#-G FILE
alternateCode=None

#	Parameter for converting between scores and likelihood ratios. This affects the column ambiguity estimates. A score is converted to a likelihood ratio by this formula: exp(score / TEMPERATURE). The default value is 1/lambda, where lambda is the scale factor of the scoring matrix, which is calculated by the method of Yu and Altschul (YK Yu et al. 2003, PNAS 100(26):15688-93).
#-t TEMPERATURE
scoreToLikelihood=None

#This option affects gamma-centroid and LAMA alignment only.
#Gamma-centroid alignments minimize the ambiguity of paired letters. In fact, this method aligns letters whose column error probability is less than GAMMA/(GAMMA+1). When GAMMA is low, it aligns confidently-paired letters only, so there tend to be many unaligned letters. When GAMMA is high, it aligns letters more liberally.
#LAMA (Local Alignment Metric Accuracy) alignments minimize the ambiguity of columns (both paired letters and gap columns). When GAMMA is low, this method produces shorter alignments with more-confident columns, and when GAMMA is high it produces longer alignments including less-confident columns.
#In summary: to get the most accurately paired letters, use gamma-centroid. To get accurately placed gaps, use LAMA.
#Note that the reported alignment score is that of the ordinary gapped alignment before realigning with gamma-centroid or LAMA.
#-g GAMMA
gammaVal=None

#Output type: 0 means counts of initial matches (of all lengths); 1 means gapless alignments; 2 means gapped alignments before non-redundantization; 3 means gapped alignments after non-redundantization; 4 means alignments with ambiguity estimates; 5 means gamma-centroid alignments; 6 means LAMA alignments; 7 means alignments with expected counts.
#If you use -j0, lastal will count the number of initial matches, per length, per query sequence. Options -l and -L will set the minimum and maximum lengths, and -m will be ignored. If you compare a large sequence to itself with -j0, it's wise to set option -L.
#If you use j>3, each alignment will get a "fullScore" (also known as "forward score" or "sum-of-paths score"). This is like the score, but it takes into account alternative alignments.
#If you use -j7, lastal will print an extra MAF line starting with "c" for each alignment. The first 16 numbers on this line are the expected counts of matches and mismatches: first the count of reference As aligned to query As, then the count of reference As aligned to query Cs, and so on. For proteins there will be 400 such numbers. The final ten numbers are expected counts related to gaps. They are:
#	The count of matches plus mismatches. (This may exceed the total of the preceding numbers, if the sequences have non-ACGT letters.)
#	The count of deleted letters.
#	The count of inserted letters.
#	The count of delete opens (= count of delete closes).
#	The count of insert opens (= count of insert closes).
#	The count of adjacent pairs of insertions and deletions.
#The final four numbers are always zero, unless you use generalized affine gap costs. They are:
#	The count of unaligned letter pairs.
#	The count of unaligned letter pair opens (= count of closes).
#	The count of adjacent pairs of deletions and unaligned letter pairs.
#	The count of adjacent pairs of insertions and unaligned letter pairs.
#-j NUMBER
output=None

#This option allows lastal to use sequence quality scores, or PSSMs, for the queries. 0 means read queries in fasta format (without quality scores); 1 means fastq-sanger format; 2 means fastq-solexa format; 3 means fastq-illumina format; 4 means prb format; 5 means read PSSMs. (Warning: Illumina data is not necessarily in fastq-illumina format; it is often in fastq-sanger format.)
#Warning: lastal cannot directly calculate E-values for PSSMs. The E-values (and the default value of -y) are determined by the otherwise-unused match and mismatch scores (options -r -q and -p). There is evidence these E-values will be accurate if the PSSM is "constructed to the same scale" as the match/mismatch scores (SF Altschul et al. 1997, NAR 25(17):3389-402).
#-Q NUMBER
qualityTreatment=None
