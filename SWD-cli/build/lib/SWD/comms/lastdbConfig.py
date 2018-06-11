#!/usr/bin/env python


#This file allows the user to specify options for lastdb. Each option is explained via its documentation entry from http://last.cbrc.jp/doc/lastdb.html and then assigned a variable name. Simply change the value of the variable to specify an option.
#All values should be character strings, input in quotes. For flags, such as "-c", just input a space in quotes (" ")
###### MAIN OPTIONS #####

# Specify lowercase-marking of repeats, by two digits (e.g. "-R01"), with the following meanings.
# First digit:
# Convert the input sequences to uppercase while reading them.
# Keep any lowercase in the input sequences.
# Second digit:
# Do not check for simple repeats.
# Convert simple repeats (e.g. cacacacacacacacac) to lowercase. This uses tantan (http://www.cbrc.jp/tantan/), which reliably prevents non-homologous alignments, unlike other repeat finders.
# Convert simple DNA repeats to lowercase, with tantan tuned for ~80% AT-rich genomes.
# -R DIGITS
repeatHandling=None

# Soft-mask lowercase letters. This means that, when we compare these sequences to some other sequences using lastal, lowercase letters will be excluded from initial matches. This will apply to lowercase letters in both sets of sequences.
# -c
softMask=None

#Specify a seeding scheme. The -m option will then be ignored. The built-in schemes are described in last-seeds.html.
#Any other NAME is assumed to be a file name. For an example of the format, see the seed files in the data directory. You can set other lastdb options on lines starting with #lastdb, but command line options override them. You can also set lastal options on lines starting with #lastal, which are overridden by options from a scoring scheme or the lastal command line.
#-u NAME
seedScheme=None

##### ADVANCED OPTIONS #####

#Allow initial matches to start only at every STEP-th position in each of the sequences given to lastdb (positions 0, STEP, 2xSTEP, etc). This reduces the memory usage of lastdb and lastal, and it makes lastdb faster. Its effect on the speed and sensitivity of lastal is not entirely clear.
#-w STEP
stepSize=None

# Allow initial matches to start only at positions that are "minimum" in any window of SIZE consecutive positions. "Minimum" means that the sequence starting here is alphabetically earliest.
# The "alphabetical" order depends on the seed pattern. The letter order is determined by the order of the letter groups, and letters in the same group are considered equivalent.
# The fraction of positions that are "minimum" is roughly: 2 / (SIZE + 1).
# -W SIZE
window=None

# Limit memory usage, by splitting the output files into smaller "volumes" if necessary. This will limit the memory usage of both lastdb and lastal, but it will make lastal slower. It is also likely to change the exact results found by lastal.
# Each volume will have at most BYTES bytes. (Roughly. The more masked letters or DNA "N"s, the more it will undershoot.) You can use suffixes K, M, and G to specify KibiBytes, MebiBytes, and GibiBytes (e.g. "-s 5G").
# However, the output for one sequence is never split. Since the output files are several-fold bigger than the input (unless you use -w), this means that mammalian chromosomes cannot be processed using much less than 2G (unless you use -w).
# There is a hard upper limit of about 4 billion sequence letters per volume. Together with the previous point, this means that lastdb will refuse to process any single sequence longer than about 4 billion.
# -s BYTES
volumeSize=None

# Specify the input format. 0 means fasta, 1 means fastq-sanger, 2 means fastq-solexa, and 3 means fastq-illumina. The fastq formats provide sequence quality data, which will be stored by lastdb and then used by lastal. These formats are described in lastal.html.
# -Q NUMBER
inputFormat="0"

# Divide the work between this number of threads running in parallel. 0 means use as many threads as your computer claims it can handle simultaneously. Currently, multi-threading is used for tantan masking only.
# -P THREADS
numThreads=None

# Specify a spaced seed pattern, for example "-m 110101". In this example, mismatches will be allowed at every third and fifth position out of six in initial matches.
# This option does not constrain the length of initial matches. The pattern will get cyclically repeated as often as necessary to cover any length.
# Although the 0 positions allow mismatches, they exclude non-standard letters (e.g. non-ACGT for DNA). If option -c is used, they also exclude lowercase letters.
# You can also specify transition constraints, e.g "-m 100TT1". In this example, transitions (but not transversions) will be allowed at every fourth and fifth position out of six. Alternatively, you can use Iedera's notation, for example "-m '#@#--##--#-#'" ('#' for match, '@' for transition, '-' or '_' for mismatch).
# You can specify multiple patterns by separating them with commas and/or using "-m" multiple times.
# -m PATTERN
seedPattern=None

# Specify your own alphabet, e.g. "-a 0123". The default (DNA) alphabet is equivalent to "-a ACGT". The protein alphabet (-p) is equivalent to "-a ACDEFGHIKLMNPQRSTVWY". Non-alphabet letters are allowed in sequences, but by default they are excluded from initial matches and get the mismatch score when aligned to anything. If -a is specified, -p is ignored.
# -a SYMBOLS
alphabet=None

# This option makes lastdb faster, at the expense of limiting your options with lastal. If you use (say) "-i 10", then you cannot use lastal with option m < 10.
# -i MATCHES
matches=None

# Specify the depth of "buckets" used to accelerate initial match finding. Larger values increase the memory usage of lastdb and lastal, make lastal faster, and have no effect on lastal's results. The default is to use the maximum depth that consumes at most one byte per possible match start position.
# -b DEPTH
bucketDepth=None

# Specify the type of "child table" to make: 0 means none, 1 means byte-size (uses a little more memory), 2 means short-size (uses somewhat more memory), 3 means full (uses a lot more memory). Choices > 0 make lastal a bit faster, but make lastdb slower, and have no effect on lastal's results. Some tests suggest that -C2 is a good choice: faster than -C1 and no slower than -C3.
# -C NUMBER
childTable=None

# Just count sequences and letters. This is much faster. Letter counting is never case-sensitive.
# -x
letterCounting=None

# Be verbose: write messages about what lastdb is doing.
# -v
verbose=None

# Show version information, and exit.
# -V, --version
version=None
