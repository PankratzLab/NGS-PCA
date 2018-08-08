#!/usr/bin/env bash

## 1. Installing anaconda:

## See https://docs.anaconda.com/anaconda/install/linux for more details

## download linux anaconda package for python 2.7

# wget "https://repo.anaconda.com/archive/Anaconda2-5.1.0-Linux-x86_64.sh"

# bash Anaconda2-5.1.0-Linux-x86_64.sh

## 2. Configure bioconda:

# conda config --add channels defaults
# conda config --add channels conda-forge
# conda config --add channels bioconda


## 3. Install MosDepth:
## See https://github.com/brentp/mosdepth for more details

# conda install mosdepth



## 4. Generate MosDepth script:

# where output will be written
mosDepthResultsDir=$HOME/mosdepthOutput/
# create output directory
mkdir -p "$mosDepthResultsDir"
# file which will contain commands to run (this will be generated below)
cmdFile=$mosDepthResultsDir"mosdepthCommands.txt"
# TOPMed reference genome
ref=$HOME/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
# Location of TOPMed Crams
inDir=$HOME/TOPMedCrams/

# clears $cmdFile if already present
echo "" > "$cmdFile"

#If .crams are listed in a one-column file instead, could also do something like:
## for file in $(cat fileOfCrams.txt); do

for file in $inDir*.cram; do

   out=$(basename "$file" .cram)
   echo "/path/to/mosdepth -n -t 2 --by 1000 --fasta $ref $mosDepthResultsDir$out.mos $file" >> "$cmdFile"

   # /path/to/mosdepth is typically $HOME/anaconda2/bin/mosdepth
   # "-t 2" use two threads 
   # "--by 1000" compute coverage on 1000bp bins
   # "-n" don't output per-base depth (save space and time)
   # Could add "--chrom chr1" to limit to chr1 only (~17 mins per sample for whole genome, ~ 1-2 mins for chr1 )

done

## 5. Run MosDepth script:

## execute commands in $cmdFile, possibly with gnu parallel, https://www.gnu.org/software/parallel/
# parallel --jobs 12 < "$cmdFile"
## "parallel --jobs 12 < $cmdFile" as set up would use 24 threads (2 threads per sample, 12 samples at once) 


