#!/usr/bin/env python
# coding=utf8

"""
%prog <BAM> <BED> <TXT> <INT> <BAM,...> [options]

Based on the input <BAM> file DE-peaks will be simulated in the provided <BED> regions, <TXT> chromosmes to be used and their length tab separated, 
<INT> defines the number of simulated replicates of the two samples,
the input/control Bam files <BAM,...> comma-separated, will be used to create an input/control Bam file for the simulation.
@author: Thomas Eder
"""

import matplotlib.pyplot as plt
import numpy as np
import pysam
import os
import sys

from scipy.stats import truncnorm
from future.utils import lrange
from random import random
from pybedtools import BedTool
from optparse import OptionParser

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

#main function
if __name__ == '__main__':
	
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	parser.add_option("-c", "--chrom", default='chr19', dest="chromosome", type="string", help="Chromosome used for simulation [default: %default]")
	parser.add_option("-r", "--read_length", default=50, dest="read_length", type="int", help="Read length [default: %default]")
	
	parser.add_option("-b","--beta", default=[0.5, 0.5], dest="beta_values", type="string", action='callback', callback=_callback_list_float, help="Alpha and Beta of Beta-distribution [default: %default]")
	parser.add_option("--max_reads", default=1, dest="max_reads", type="float", help="Maximum percentage of reads for repcilate numbers [default: %default]")
	parser.add_option("--min_reads", default=0, dest="min_reads", type="float", help="Minimum percentage of reads for repcilate numbers [default: %default]")
	
	parser.add_option("--rep_sd", default=0.1, dest="rep_sd", type="float", help="Standard deviation of replicate read numbers [default: %default]")
	parser.add_option("--rep_mean", default=0.9, dest="rep_mean", type="float", help="Mean of replicate read numbers [default: %default]")
	
	parser.add_option("-d", "--dp-thres", default=0.7, dest="dp_thres", type="float", help="Threshold of reads to define a DB peak [default: %default]")
	parser.add_option("-m", "--min-counts", default=10, dest="min_counts", type="int", help="Minimum number of reads for a DB peak [default: %default]")
	
	parser.add_option("-n", "--non_peak_percetage", default=0.9, dest="noise_treshold", type="float", help="Percentage of reads outside the defined peak regions [default: %default]")

	parser.add_option("--control_read_num", default=0, dest="num_control_reads", type="int", help="Number of reads in the input control Bam files, if 0 use mean of input control files[default: %default]")
	
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

	#normal distributions sampling with max arround num_reads_sample to get num_reads_replicate
	get_num_reads_rep = get_truncated_normal(mean=options.rep_mean , sd=options.rep_sd, low=options.min_reads, upp=options.max_reads)
	num_reads_s1rep = list(get_num_reads_rep.rvs(num_replicates))
	get_num_reads_rep = get_truncated_normal(mean=options.rep_mean , sd=options.rep_sd, low=options.min_reads, upp=options.max_reads)
	num_reads_s2rep = list(get_num_reads_rep.rvs(num_replicates))
	
	
	#read bed
	input_bed = BedTool(inbed_filename)
	#open output bed
	output_bed = open(options.outbed_filename+".bed", "w")
	#header
	output_bed.write("#Chromsome\tStart\tStop\tSample1_mean_counts\tSample2_mean_counts%s%s\tStatus\tPeak_name\n" %( \
	"".join(list(map(lambda x: '\tSample1_replicate%i_counts' % x, lrange(1,num_replicates+1)))),\
	"".join(list(map(lambda x: '\tSample2_replicate%i_counts' % x, lrange(1,num_replicates+1))))))

	region_count=1

	#for each bed entry get reads...
	for region in input_bed:
		#counts per replicate
		sample1_rep_count = np.zeros(num_replicates)
		sample2_rep_count = np.zeros(num_replicates)
		
		#beta distribution to get aprox. same size or down regulated read number for sample x, here we deside which sample to prefer
		read_frac = np.random.beta(options.beta_values[0], options.beta_values[1])
		
		#for each of the regions change height based on beta distribution
		for read in input_bam.fetch(region.chrom, region.start, region.end):
			
			if random() <= read_frac: #add to sample 1 or sample 2 and to list
				rep_cnt=0
				#now decide to which replicate
				for rep_treshold in num_reads_s1rep:
					if random() <= rep_treshold:
						output_s1_bams[rep_cnt].write(read)
						sample1_rep_count[rep_cnt]+=1
					
					rep_cnt+=1
			else:
				rep_cnt=0
				#now decide to which replicate
				for rep_treshold in num_reads_s2rep:
					if random() <= rep_treshold:
						output_s2_bams[rep_cnt].write(read)
						sample2_rep_count[rep_cnt]+=1
					
					rep_cnt+=1

		#get mean count per sample for this region
		s1_mean_count=0
		s1_count=0
		s2_mean_count=0
		s2_count=0
		
		for count in sample1_rep_count:
			s1_count+=count
		s1_mean_count=s1_count/num_replicates
		
		for count in sample2_rep_count:
			s2_count+=count
		s2_mean_count=s2_count/num_replicates
		
		#is it DE?
		stat=0
		if s1_mean_count > s2_mean_count and s1_mean_count > options.min_counts and (s1_mean_count / (s1_mean_count + s2_mean_count)) >= options.dp_thres: #up in sample1
			stat=1
		elif s2_mean_count > s1_mean_count and s2_mean_count > options.min_counts and (s2_mean_count / (s1_mean_count + s2_mean_count)) >= options.dp_thres: #up in sample2
			stat=2
		else: #not meeting user defined DP parameters or is not up in either sample
			stat=0
			
		#id
		name='peak'+str(region_count)
		
		#now save data for this region
		output_bed.write("%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n" %(region.chrom, region.start, region.end, s1_mean_count, s2_mean_count,\
		"\t".join(list(map(lambda x: str(x), sample1_rep_count))),  "\t".join(list(map(lambda x: str(x), sample2_rep_count))), stat, name))
		
		region_count+=1


	#now take all other regions and leave them as they are...
	complement_bed = input_bed.complement(g=in_chromsizes)
	
	for region in complement_bed:
		
		#check region if its to narrow
		if not (region.start+options.read_length+(2+options.read_length) >= region.end-options.read_length): #leave region out if to narrow we alread processes these reads...
			
			for read in input_bam.fetch(region.chrom, region.start+options.read_length, region.end-options.read_length): #decrease size of region by read length as reads with only 1 bp match will also be fetched for a spezific region
				
				if random() <= options.noise_treshold: #do not just take all, bring in some variance
					rep_cnt=0
					for rep_treshold in num_reads_s1rep:
						if random() <= rep_treshold:
							output_s1_bams[rep_cnt].write(read)
						rep_cnt+=1
					
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
		options.num_control_reads = len(all_control_reads)/len(input_control_bam_files)

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
