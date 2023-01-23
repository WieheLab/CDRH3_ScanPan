# CDRH3_ScanPan

The python code for generating the distance counts from a sequence set generated either experimentally or *in silico*. The CDRH3_ScanPan.py program is designed to take in sets of sequences and scoring matrix. We use a scoring matrix generated from single amino acid substitutions using mutagenesis libraries.

CDRH3_python_requirements.txt is a list of nonstandard python libraries needed to run CDRH3_ScanPan.py.

## Running CDRH3_ScanPan.py

for the help menu
```
python CDRH3_ScanPan.py -h
```

to run:
```
python CDRH3_ScanPan.py -m <matrixfile> -s <seqset>
```

## Arguments

-m \<matrix\>  -  scoring matrix to use for scoring sequences

-s \<sequence file\>  -  file containing all the sequences, expecting fasta format

-d \<distance file\>  -  name of generated distance file

-g \<gene\>  -  looks for specific D gene

-c \<cutoff\>  -  cut off used to determine a hit or miss from the matrix file when calculating the distance

-o \<pdffile name\>  -  output pdf name for the results

--scale <scale to use>  -  scale to use when generating the histograms with the default set to the number of sequences

--writetop <write seqs with <5 dist to file>

--noseqout <do not write seqs with <5 dist to file>


# Directory Generate_Igor_Seqs

The directory Generate_Igor_Seqs contains scripts for generating random sets of sequences using the IGoR software (https://statbiophys.github.io/IGoR/). The IGoR software is required to run these scripts. This is broken into job scripts (for Slurm) and bash scripts written to be run on the Duke computing cluster.

**igor_runs.job** - is the slurm script used to generate random heavy sequences using IGoR\
**igor_run_forCounts.job** - A modification of the igor_runs.job to return number of sequences generated\
**collect_sort.job** - takes the generated sequences and uses weblogo to create logo plots\
**collect.job** - Job script to take the generated sequences from igor_runs.job (or igor_run_forCounts.job) and put them into one file and calls create_cdr3_seqs.sh if not enough unique sequences are generated\
**create_cdr3_seqs.sh** - simple bash script to call igor_runs.job and when completed the collect.job script


# Directory Scripts

run_weblogo.V3.sh - bash script to generate logo plots using weblogo \
calc_SE.job - bash script job script for running the Calc_Site_Variation.py script on multiple sequences\
Calc_Site_Variation.py - python script used to calculate the shannon entropy of a sequence set\
compress.job - bash script job script for slurm to use gzip to compress consensus.txt files\
