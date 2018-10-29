#!/bin/bash

## Generate mappability track using https://github.com/gemtools/gemtools

# pref=$root 
# reference=$refGenome 
# idxpref="${pref}_index"
# thr=10; #change it to the number you want to use
# outmappa=$root".mappa.tsv" # should be changed
# #make index for the genome
# gem-indexer -T ${thr} -c dna -i ${reference} -o ${idxpref}

# # The following calculates index and creates mappability tracks with various kmer lengths.
# # it may take long time to finish.
# # choose the right kmer length for you.
# for kmer in 100 125 150; do

#   # compute mappability data
#   gem-mappability -T ${thr} -I ${idxpref}.gem -l ${kmer} -o ${pref}_${kmer}
#   mpc=$(gawk '{c+=gsub(s,s)}END{print c}' s='!' ${pref}_${kmer}.mappability)
#   echo ${pref}_${kmer}"\t"$mpc >> $outmappa
#   # convert results to wig and bigwig
#   gem-2-wig -I ${idxpref}.gem -i ${pref}_${kmer}.mappability -o ${pref}_${kmer}
#   /home/pankrat2/shared/bin/hgDownload/wigToBigWig ${pref}_${kmer}.wig ${pref}.sizes ${pref}_${kmer}.bw
#   wig2bed < ${pref}_${kmer}.wig > ${pref}_${kmer}.bed
# done

## Extract regions of poor mappability


threshold=0.75

# Pull out regions with mappability < threshold
awk -v threshold="$threshold" '$5 < threshold {print}' GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_150.bed |cut -f 1-3 > GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_150.$threshold.bed

bedtools merge -i GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_150.$threshold.bed >GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_150.$threshold.merge.bed


## Download sv_blacklist.bed http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed

cut -f 1-3 sv_blacklist.bed >ngs_pca_exclude.bed
cat GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_150.$threshold.merge.bed >>ngs_pca_exclude.bed

sort -k1,1 -k2,2n ngs_pca_exclude.bed > ngs_pca_exclude.sorted.bed 

bedtools merge -i ngs_pca_exclude.sorted.bed > ngs_pca_exclude.sorted.merge.bed 

gzip ngs_pca_exclude.sorted.merge.bed 