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

threshold=1.0

 # 25 100 125 150
for kmer in 35 50; do


# Pull out regions with mappability < threshold
awk -v threshold="$threshold" '$5 < threshold {print}' GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_$kmer.bed |cut -f 1-3 > GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_$kmer.$threshold.bed

bedtools merge -i GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_$kmer.$threshold.bed >GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_$kmer.$threshold.merge.bed


## Download sv_blacklist.bed http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed

cut -f 1-3 sv_blacklist.bed >ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.bed

cat GRCh38_full_analysis_set_plus_decoy_hla.chr1-chr22-X-Y-M_$kmer.$threshold.merge.bed >>ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.bed

sort -k1,1 -k2,2n ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.bed > ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.sorted.bed

bedtools merge -i ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.sorted.bed >  ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.sorted.merge.bed

# cp ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.sorted.merge.bed ~/git/NGS-PCA/resources

# gzip ~/git/NGS-PCA/resources/ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.sorted.merge.bed


#Download DGV database wget "http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2016-08-31.txt"

awk '{print "chr"$2"\t"$3"\t"$4}' GRCh38_hg38_variants_2016-08-31.txt |grep -v "start" > ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.bed
cat ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.sorted.merge.bed >> ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.bed

sort -k1,1 -k2,2n ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.bed > ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.bed

bedtools merge -i ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.bed > ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.merge.bed

cp ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.merge.bed ~/git/NGS-PCA/resources

gzip ~/git/NGS-PCA/resources/ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.merge.bed

cat ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.merge.bed > ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.bed
cut -f 1-3 /Users/Kitty/git/Analysis/mappability/genomicSuperDups.hg38.bed >> ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.bed

sort -k1,1 -k2,2n ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.bed | bedtools merge > ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed
cp ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed ~/git/NGS-PCA/resources
gzip ~/git/NGS-PCA/resources/ngs_pca_exclude.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed


done

