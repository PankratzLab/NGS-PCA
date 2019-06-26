#!/bin/bash

threshold=1
kmer=50
# Pull out regions with mappability < threshold
awk -v threshold="$threshold" '$5 < threshold {print}' hg19_canonical_50.bed |cut -f 1-3 |bedtools merge |bgzip > hg19_canonical_$kmer.$threshold.merge.bed.gz
tabix hg19_canonical_$kmer.$threshold.merge.bed.gz


## Download sv_blacklist.bed wget "http://cf.10xgenomics.com/supp/genome/hg19/sv_blacklist.bed"

cut -f 1-3 ~/git/NGS-PCA/resources/hg19/sv_blacklist.bed >ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.bed

cat hg19_canonical_$kmer.$threshold.merge.bed.gz |zcat  >> ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.bed


bedtools sort -i ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.bed|bedtools merge > ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.sorted.merge.bed

#Download DGV database wget "http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt"

awk '{print "chr"$2"\t"$3"\t"$4}' GRCh37_hg19_variants_2016-05-15.txt |grep -v "start" > ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.bed
cat ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.sorted.merge.bed >> ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.bed

bedtools sort -i ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.bed | bedtools merge > ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.merge.bed


#Download genomic super dups
# wget "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=732005553_B6GUBYxSG2pHft9qHx0AzBfdilE6&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_genomicSuperDups&hgta_ctDesc=table+browser+query+on+genomicSuperDups&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED"
# mv "hgTables?hgsid=732005553_B6GUBYxSG2pHft9qHx0AzBfdilE6&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_genomicSuperDups&hgta_ctDesc=table+browser+query+on+genomicSuperDups&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbD" genomicSuperDups.hg19.bed

cat ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.sorted.merge.bed > ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.bed
cut -f 1-3 genomicSuperDups.hg19.bed >> ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.bed

bedtools sort -i ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.bed | bedtools merge |bgzip > ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed.gz
cp ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed.gz /Users/Kitty/git/NGS-PCA/resources/hg19/ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed.gz
tabix /Users/Kitty/git/NGS-PCA/resources/hg19/ngs_pca_exclude.hg19.sv_blacklist.map.kmer.$kmer.$threshold.dgv.gsd.sorted.merge.bed.gz



