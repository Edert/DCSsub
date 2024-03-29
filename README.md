# Sub-sampling of differential ChIP-seq data

DCSsub sub-samples reads from one or more genuine ChIP-seq BAM files for two samples (e.g. treatment and control) in two or more replicates. We strongly encourage new users to start with test-fragement-count-distribution.py. This script should be used to test the parameter settings without creating fasta output files but creating informative plots.

## DCSsub ##
**DCSsub.py** \<BAM\> \<BED\> \<TXT\> \<INT\> \<BAM,...\> [options]

Based on the input \<BAM\> file(s), differential-peak regions will be sub-sampled in the provided \<BED\> regions, a \<TXT\> file with a list of chromosomes to be used and their length tab separated, \<INT\> defines the number of simulated replicates of the two samples, the input/control Bam files \<BAM,...\> comma-separated, will be used to create an input/control Bam file for the simulation.

Options:

-c, --chrom	Chromosome used for sub-sampling, default='chr19'

-r, --read_length	Read length, default=50

-b, --beta	Alpha and Beta of Beta-distribution, default=[0.5, 0.5]

--max_reads	[Non peak region] Maximum percentage of reads for replicates, default=1

--min_reads	[Non peak region] Minimum percentage of reads for replicates, default=0

--rep_sd	[Non peak region] Standard deviation of replicate read numbers, default=0.1

--rep_mean	[Non peak region] Mean of replicate read numbers, default=0.9

--frag-count-scaling	Scaling of frag distribution, no scaling, scaling of beta result based on read counts (with exponential distribution) or scaling of read counts based on beta result (with Laplace distribution): none, frag, beta, default="none"

--frag-count-lp-scale	Scale for Laplace distribution if frag-count-scaling is frag, default=0.1

--frag-count-ex-loc	Loc for exponential distribution if frag-count-scaling is beta, default=10

--frag-count-ex-scale	Scale for exponential distribution if frag-count-scaling is beta, default=100

-d, --dp-thres	Threshold of reads to define a DB peak , default=0.7

-m, --min-counts	Minimum number of reads for a DB peak, default=10

-n, --non_peak_percetage	Percentage of reads outside the defined peak regions , default=0.9

-s, --skewness	Variance between replicates (the higher, the less variance), default=10

--control_read_num	Number of reads in the input control Bam files, if 0 use mean of input control files, default=0

--out_bed	Name of the results bed file, default="results.bed"

--out_sample	Basename of the results Bam files, default="out"

--out_control	Basename of the input/control Bam files, default="out"

## test-fragement-count-distribution ##

**test-fragement-count-distribution.py** \<BAM\> \<BED\> \<TXT\> \<INT\> [options]
  
Based on the input \<BAM\> file DE-peaks will be sub-sampled in the provided \<BED\> regions, \<TXT\> chromosomes to be used and their length tab separated, \<INT\> defines the number of simulated replicates of the two samples.
  
Options:
  
-c, --chrom	Chromosome used for simulation, default='chr19'
  
--frag-count-scaling	Scaling of read distribution, no scaling, scaling of beta result based on frag counts (with exponential distribution) or scaling of frag counts based on beta result (with Laplace distribution): none, frag, beta, default="none"
  
--frag-count-lp-scale	Scale for Laplace distribution if frag-count-scaling is frag, default=0.1
  
--frag-count-ex-loc	Loc for exponential distribution if frag-count-scaling is beta, default=10
  
--frag-count-ex-scale	Scale for exponential distribution if frag-count-scaling is beta, default=100
  
--beta	Alpha and Beta of Beta-distribution, default=[0.5, 0.5]
  
 
## Citation ##

If you are using this tool or parts of it in your research, please cite:

Eder, T., Grebien, F. Comprehensive assessment of differential ChIP-seq tools guides optimal algorithm selection. Genome Biol 23, 119 (2022). https://doi.org/10.1186/s13059-022-02686-y
