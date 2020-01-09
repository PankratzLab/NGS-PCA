#!/bin/bash


# wget  -nd  biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/xgen_plus_spikein.b38.bed
# http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=3801

# bgzip xgen_plus_spikein.b38.bed
# tabix xgen_plus_spikein.b38.bed.gz


cat xgen_plus_spikein.b38.bed.gz |zcat |awk '{if($2-20000 <0) print "chr"$1"\t"0"\t"$3+20000; else print "chr"$1"\t"$2-20000"\t"$3+20000}' |bgzip >xgen_plus_spikein.b38.chr.bed.gz
cat ../ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz xgen_plus_spikein.b38.chr.bed.gz|zcat |bedtools sort |bedtools merge |bgzip >ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.xgen.sorted.merge.bed.gz
tabix ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.xgen.sorted.merge.bed.gz






