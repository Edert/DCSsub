#!/usr/bin/env python
# coding=utf8

"""
%prog <BAM> <BED> <TXT> <INT> <BAM,...> [options]

Based on the input <BAM> file DE-peaks will be sub-sampled in the provided <BED> regions, <TXT> chromosmes to be used and their length tab separated, 
<INT> defines the number of simulated replicates of the two samples,
the input/control Bam files <BAM,...> comma-separated, will be used to create an input/control Bam file for the simulation.
@author: Thomas Eder
"""
from __future__ import print_function
from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pysam
import os
import sys

from scipy.stats import truncnorm
from future.utils import lrange
from random import random
from random import randrange
from pybedtools import BedTool
from optparse import OptionParser
from numpy.random import dirichlet

class HelpfulOptionParser(OptionParser):
	"""An OptionParser that prints full help on errors."""
	def error(self, msg):
		self.print_help(sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def sort_index(file_name):
	"""Sort and index Bam-file"""
	pysam.sort("-o",'tmp.bam',file_name)
	os.rename('tmp.bam',file_name)
	pysam.index(file_name)

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
	"""Get truncnorm object with defined mean and sd, it represents a truncated normal distribution"""
	return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
	
def _callback_list_float(option, opt, value, parser): #get a list of float options and save it
	try:
		tmplist = list(map(lambda x: float(x), value.split(',')))
		setattr(parser.values, option.dest, tmplist)  
	except Exception as e:
		print("got exception: %r, please check input parameters" % (e,))
		exit()

def _scale_laplace(nfrags, read_frac, scale):
	"""Scale beta result based on Laplace distribution"""
	loc = 0.5 #fixed
	
	top_nfrgs_scale = laplace.pdf(0.5, loc, scale)#make sure at postion 0.5 is 100%
	scale_factor = 1 / top_nfrgs_scale
	nfrags_scale = laplace.pdf(read_frac, loc, scale) *scale_factor #the bigger the difference from beta the smaller the sample...
	
	if(nfrags_scale == 0):
		nfrags = nfrags 
	else:
		nfrags = int(nfrags *nfrags_scale) #new fragment scaling based on beta result
		
	return nfrags

def _scale_lognorm(nfrags, read_frac, sigma, loc, scale):
	"""Scale number_frags result based on Lognorm distribution"""
	#loc=scale*-1
	
	top_beta_scale = sp.stats.lognorm.pdf(1,s=sigma,loc=loc, scale=scale)#to get to 100% at size 10
	scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at 1
	beta_scale = sp.stats.lognorm.pdf(nfrags,s=sigma,loc=loc, scale=scale) *scale_factor
	
	if(beta_scale == 0):
		read_frac = read_frac 
	else:
		read_frac = ((read_frac-0.5)*beta_scale)+0.5 #center and scale it based on beta_scale from exponential distri, so reduce to 0.5 for high values
	return read_frac

def _scale_exp(nfrags, read_frac, loc, scale):
	"""Scale number_frags result based on Exponential distribution"""
	
	top_beta_scale = sp.stats.expon.pdf(loc,loc=loc, scale=scale)#to get to 100% at size loc
	scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at 1
	beta_scale = sp.stats.expon.pdf(nfrags,loc=loc, scale=scale) *scale_factor
	
	read_frac_old = read_frac
	
	if(beta_scale == 0):
		read_frac = read_frac 
	else:
		read_frac = ((read_frac-0.5)*beta_scale)+0.5 #center and scale it based on beta_scale from exponential distribution, we do this to reduce the read_frac towards 0.5 for high values
	
	return read_frac


def _reads_to_replicates(read_list, skewness, n_replicates):
	"""puts reads in a randomly selected replicate"""
	#get reads into replicates
	bins = np.cumsum([0] + list(dirichlet(np.array([skewness] * n_replicates), 1)[0])) #sample dirichlet dist, get array, make list, add 0 to beginning and sum up all values to 1, so we end up with an array of bins with probabilities
	replicate_counts = [0] * n_replicates
	replicate_list = [[] for i in lrange(n_replicates)]
	
	for readobj in read_list: #go through reads
		index = np.digitize([random()], bins)[0].item() - 1 #get index-1 of random number (between 0 or 1) fitting to the right bin
		replicate_list[index].append(readobj) #add element to replicate via choosen index
		replicate_counts[index] += 1
	
	return replicate_list, replicate_counts #return number of reads per replicate

#main function
if __name__ == '__main__':
	
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	parser.add_option("-c", "--chrom", default='chr19', dest="chromosome", type="string", help="Chromosome used for sub-sampling [default: %default]")
	parser.add_option("-r", "--read_length", default=50, dest="read_length", type="int", help="Read length [default: %default]")
	
	parser.add_option("-b","--beta", default=[0.5, 0.5], dest="beta_values", type="string", action='callback', callback=_callback_list_float, help="Alpha and Beta of Beta-distribution [default: %default]")
	
	parser.add_option("--max_reads", default=1, dest="max_reads", type="float", help="[Non peak region] Maximum percentage of reads for replicates [default: %default]")
	parser.add_option("--min_reads", default=0, dest="min_reads", type="float", help="[Non peak region] Minimum percentage of reads for replicates [default: %default]")
	parser.add_option("--rep_sd", default=0.1, dest="rep_sd", type="float", help="[Non peak region] Standard deviation of replicate read numbers [default: %default]")
	parser.add_option("--rep_mean", default=0.9, dest="rep_mean", type="float", help="[Non peak region] Mean of replicate read numbers [default: %default]")
	
	parser.add_option("--read-count-scaling", default="none", dest="frag_count_scaling", type="string", help="Scaling of read distribution, no scaling, scaling of beta result based on read counts (with exponential distribution) or scaling of read counts based on beta result (with Laplace distribution): none, frag, beta [default: %default]")
	parser.add_option("--read-count-lp-scale", default=0.1, dest="frag_count_lp_sc", type="float", help="Scale for Laplace distribution if read-count-scaling is frag [default: %default]")
	parser.add_option("--read-count-ex-loc", default=10, dest="frag_count_ex_lo", type="float", help="Loc for exponential distribution if read-count-scaling is beta [default: %default]")
	parser.add_option("--read-count-ex-scale", default=100, dest="frag_count_ex_sc", type="float", help="Scale for exponential distribution if read-count-scaling is beta [default: %default]")
	
	parser.add_option("-d", "--dp-thres", default=0.7, dest="dp_thres", type="float", help="Threshold of reads to define a DB peak [default: %default]")
	parser.add_option("-m", "--min-counts", default=10, dest="min_counts", type="int", help="Minimum number of reads for a DB peak [default: %default]")
	
	parser.add_option("-n", "--non_peak_percetage", default=0.9, dest="noise_treshold", type="float", help="Percentage of reads outside the defined peak regions [default: %default]")
	
	parser.add_option("-s", "--skewness", default=10, dest="skewness", type="float", help="Variance between replicates (the higher, the less variance) [default: %default]")
	
	parser.add_option("--control_read_num", default=0, dest="num_control_reads", type="int", help="Number of reads in the input control Bam files, if 0 use mean of input control files [default: %default]")
	
	parser.add_option("--out_bed", default="results.bed", dest="outbed_filename", type="string", help="Name of the results bed file [default: %default]")
	parser.add_option("--out_sample", default="out", dest="output_name", type="string", help="Basename of the results Bam files [default: %default]")
	parser.add_option("--out_control", default="out", dest="output_control_name", type="string", help="Basename of the input/control Bam files [default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 5
	
	#read in required input arguments
	if len(args) != num_args:
		parser.error("At minimum %s parameters are requred to run this tool" % num_args)
	input_bam_file = args[0]
	inbed_filename = args[1]
	in_chromsizes = args[2]
	num_replicates = int(args[3])
	input_control_bam_files_string = args[4]
	input_control_bam_files = input_control_bam_files_string.split(",")

	#get sample bam file...
	input_bam = pysam.AlignmentFile(input_bam_file,'rb')
	subsample_chrom = options.chromosome 

	#set up output bams
	output_s1_bams=[]
	output_s2_bams=[]
	
	for x in lrange(1,num_replicates+1):
		output_s1_bams.append(pysam.AlignmentFile(options.output_name+"_sample1-rep"+str(x)+".bam", "wb", template=input_bam))
		output_s2_bams.append(pysam.AlignmentFile(options.output_name+"_sample2-rep"+str(x)+".bam", "wb", template=input_bam))
	
	#read bed
	input_bed = BedTool(inbed_filename)
	#open output bed
	output_bed = open(options.outbed_filename+".bed", "w")
	#header
	output_bed.write("#Chromsome\tStart\tStop\tSample1_mean_counts\tSample2_mean_counts%s%s\tStatus\tPeak_name\n" %( \
	"".join(list(map(lambda x: '\tSample1_replicate%i_counts' % x, lrange(1,num_replicates+1)))),\
	"".join(list(map(lambda x: '\tSample2_replicate%i_counts' % x, lrange(1,num_replicates+1))))))

	region_count=1
	
	#normal distributions sampling with max arround num_reads_sample to get num_reads_replicate
	get_num_reads_rep = get_truncated_normal(mean=options.rep_mean , sd=options.rep_sd, low=options.min_reads, upp=options.max_reads)
	num_reads_s1rep = list(get_num_reads_rep.rvs(num_replicates))
	get_num_reads_rep = get_truncated_normal(mean=options.rep_mean , sd=options.rep_sd, low=options.min_reads, upp=options.max_reads)
	num_reads_s2rep = list(get_num_reads_rep.rvs(num_replicates))

	#for each bed entry get reads...
	sample1_replicate_counts = [0] * num_replicates
	sample2_replicate_counts = [0] * num_replicates
	sample1_replicate_list = [[] for i in lrange(num_replicates)]
	sample2_replicate_list = [[] for i in lrange(num_replicates)]
	multi_factor = 2*num_replicates #per read multiply its use, so that we can definitely fill 2 samples with x replicates
	
	for region in input_bed:
		sample1_list = []
		sample2_list = []
		sample1_count = 0
		sample2_count = 0
		
		#beta distribution to get aprox. same size or down regulated read number for sample x, here we deside which sample to prefer
		read_frac = np.random.beta(options.beta_values[0], options.beta_values[1])
		
		#get number of reads in this region, name id nfragments like in DCSsim
		number_frags = input_bam.count(region.chrom, region.start, region.end)
		
		#now choose which way we want to select the number of fragments
		if( options.frag_count_scaling == "none"):
			#no scaling, no changes in nfrags and beta
			number_frags = number_frags
			read_frac = read_frac
		elif (options.frag_count_scaling == "beta") :#nfrags scaling via Laplace based on beta result, beta not changed
			number_frags = _scale_laplace(number_frags, read_frac, options.frag_count_lp_sc)
		elif (options.frag_count_scaling == "frag") :#nfrags scaling via exp distribution based on fragment counts, number_frags not changed
			read_frac = _scale_exp(number_frags, read_frac, options.frag_count_ex_lo, options.frag_count_ex_sc)
		else:
			print("Unknown scaling method, %s, please choose 'none','frag' or 'beta', exiting now" % (frag_count_scaling))
			exit(1)
		
		#for each region change height based on beta distribution
		read_counter=0
		read_pool=[]
		
		#first get all reads of this region
		for read in input_bam.fetch(region.chrom, region.start, region.end):
			read_pool.append(read)
		
		#then randomly select on of the reads for x-times
		for _ in lrange(number_frags):
			for _ in lrange(multi_factor): #do this x times
				if(read_counter > number_frags): #check if we meet calculated and scaled (if turned on) number of fragments/reads
					break
				read_counter+=1
				
				#get one random read out of region-pool
				random_index = randrange(len(read_pool))
				read = read_pool[random_index]
				
				if random() <= read_frac: #add to sample 1 or sample 2 and to list
					sample1_count += 1
					sample1_list.append(read)
				else:
					sample2_count += 1
					sample2_list.append(read)
			else:
				continue
			break
		
		#now decide to which replicate
		sample1_replicate_list, sample1_replicate_counts = _reads_to_replicates(sample1_list, options.skewness, num_replicates)
		sample2_replicate_list, sample2_replicate_counts = _reads_to_replicates(sample2_list, options.skewness, num_replicates)
		
		#write reads to files
		for repnum in lrange(num_replicates):
			for read in sample1_replicate_list[repnum]:
				output_s1_bams[repnum].write(read)
			for read in sample2_replicate_list[repnum]:
				output_s2_bams[repnum].write(read)
		
		#get mean
		sample1_mean_count = sample1_count / num_replicates
		sample2_mean_count = sample2_count / num_replicates
		
		#check if me meet threshold
		if sample1_mean_count > sample2_mean_count and sample1_mean_count > options.min_counts and (sample1_mean_count / (sample1_mean_count + sample2_mean_count)) >= options.dp_thres: #up in sample1
			stat=1
		elif sample2_mean_count > sample1_mean_count and sample2_mean_count > options.min_counts and (sample2_mean_count / (sample1_mean_count + sample2_mean_count)) >= options.dp_thres: #up in sample2
			stat=2
		else: #not meeting user defined DP parameters or is not up in either sample
			stat=0
		
		#id
		name='peak'+str(region_count)
		
		#now save data for this region
		output_bed.write("%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n" %(region.chrom, region.start, region.end, sample1_mean_count, sample2_mean_count,\
		"\t".join(list(map(lambda x: str(x), sample1_replicate_counts))),  "\t".join(list(map(lambda x: str(x), sample2_replicate_counts))), stat, name))
		
		region_count+=1
	
	
	#now take all other regions and leave them as they are...
	complement_bed = input_bed.complement(g=in_chromsizes)
	
	for region in complement_bed:
		
		#check region if its to narrow
		if not (region.start+options.read_length+(2+options.read_length) >= region.end-options.read_length): #leave region out if to narrow we alread processes these reads...
			
			for read in input_bam.fetch(region.chrom, region.start+options.read_length, region.end-options.read_length): #decrease size of region by read length as reads with only 1 bp match will also be fetched for a spezific region
				#sample1
				if random() <= options.noise_treshold: #do not just take all, bring in some variance
					rep_cnt=0
					for rep_treshold in num_reads_s1rep:
						if random() <= rep_treshold:
							output_s1_bams[rep_cnt].write(read)
						rep_cnt+=1
				#sample2
				if random() <= options.noise_treshold:
					rep_cnt=0
					for rep_treshold in num_reads_s2rep:
						if random() <= rep_treshold:
							output_s2_bams[rep_cnt].write(read)
						rep_cnt+=1


	#input/control bam file(s)
	all_control_reads=[]

	#read all and store
	for input_control_bam_file in input_control_bam_files:
		input_control_bam = pysam.AlignmentFile(input_control_bam_file,'rb')
		for read in input_control_bam.fetch(options.chromosome):
			all_control_reads.append(read)
		input_control_bam.close()

	if options.num_control_reads == 0 : #if not spezified by user use mean number of reads in the input control files
		options.num_control_reads = len(all_control_reads)//len(input_control_bam_files)

	#select random and write one control file per sample 
	selected_control_reads = np.random.choice(all_control_reads, options.num_control_reads, replace=False)
	output_s1_control_bam = pysam.AlignmentFile(options.output_control_name+"_sample1-INPUT.bam", "wb", template=input_bam)
	for read in selected_control_reads:
		output_s1_control_bam.write(read)

	selected_control_reads = np.random.choice(all_control_reads, options.num_control_reads, replace=False)
	output_s2_control_bam = pysam.AlignmentFile(options.output_control_name+"_sample2-INPUT.bam", "wb", template=input_bam)
	for read in selected_control_reads:
		output_s2_control_bam.write(read)


	#close bam files and index new bam files
	input_bam.close()
	output_bed.close()

	for out_bam in output_s1_bams:
		out_bam.close()
		sort_index(out_bam.filename)

	for out_bam in output_s2_bams:
		out_bam.close()
		sort_index(out_bam.filename)

	#close, sort and index control files
	output_s1_control_bam.close()
	sort_index(output_s1_control_bam.filename)

	output_s2_control_bam.close()
	sort_index(output_s2_control_bam.filename)
	
	exit(0)
