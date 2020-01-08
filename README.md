# NGS-PCA
Methods for running PCA on NGS data


## Install [mosdepth](https://github.com/brentp/mosdepth)

Install mosdepth using the instructions from https://github.com/brentp/mosdepth#installation.
There are lots of ways to do this, including downloading a linux executable, as a Docker image, or using brew.

## Run mosdepth on bams or crams
The "--by 1000" (compute coverage on 1000bp bins) is really the only important argument, and each run is going to look something like this:
`mosdepth -n -t 1 --by 1000 --fasta /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa output_filename input_filename.bam`

but here is an example script to iterate over all BAM files in a directory (could be run on CRAM files as well

```bash
ref=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa

dirOfBams=/path/to/bams/
mosdepthResultsDir=/path/to/mosdepthOutput/
mosdepthThreads=1
parallelThreads=24

find "$dirOfBams" -type f -name "*.bam" \
|parallel -j $parallelThreads "mosdepth -n -t $mosdepthThreads --by 1000 --fasta $ref $mosdepthResultsDir{/.}.by1000 {}"
```

## Run ngs pca

The number of PCs to compute should be in the range of 5% of your sample size and likely far more than you'll actually use.
We're still working on optimizing the number of iterations and oversample parameter - but this should be reasonable.
This will generate svd.pcs.txt in the output directory


```bash
ngsPCAOutputDir=/path/to/ngsPCA/
ngsPCAThreads=24
# number of PCs to compute, will likely only use ~10 for 700 samples, so computing 100 should be plenty to play with
numPCs=100

ngsPCAExcludeRegions=ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed
jar=$HOME/ngspca-0.01-SNAPSHOT.846fb69.jar

java -Xmx60G -jar "$jar" \
-input $mosDepthResultsDir \
-outputDir $ngsPCAOutputDir \
-numPC $numPCs \
-sampleEvery 0 \
-threads $ngsPCAThreads \
-iters 	10 \
-randomSeed 42 \
-oversample 100 \
-bedExclude $ngsPCAExcludeRegions

```
`ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed` can be found [here](https://github.com/PankratzLab/NGS-PCA/blob/master/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.50.1.0.dgv.gsd.sorted.merge.bed.gz)

The jar can be downloaded from the the latest release https://github.com/PankratzLab/NGS-PCA/releases or https://github.com/PankratzLab/NGS-PCA/blob/master/build/ngspca-0.01-SNAPSHOT.846fb69.jar

The ngspca jar will essentially:

1. Select autosomal bins that do not overlap any region in the excluded bed
2. Normalize input data ( https://github.com/PankratzLab/NGS-PCA )
	- Normalize within sample by computing fold change 
		- log2(coverage of bin / median coverage of all selected bins)
	- Center each bin to median fold-change of 0 across all samples 
3. Perform Randomized PCA
	- Described in https://epubs.siam.org/doi/abs/10.1137/090771806 and https://epubs.siam.org/doi/abs/10.1137/100804139
  	- Similar to the https://github.com/erichson/rSVD R package

