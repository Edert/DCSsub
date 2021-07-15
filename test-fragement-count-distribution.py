#!/usr/bin/env python
# coding=utf8

"""
%prog <BAM> <BED> <TXT> <INT> [options]

Based on the input <BAM> file DE-peaks will be sub-sampled in the provided <BED> regions, <TXT> chromosmes to be used and their length tab separated, 
<INT> defines the number of simulated replicates of the two samples

Creates a histogram of all fragments count distribution based on regions from <BED> of <INT> replicates for two samples, for all samples, for sample1 and sample2 and creates an MA-plot based on the supplied beta-values

@author:  Thomas Eder
"""

from __future__ import print_function
from __future__ import division

import sys
import math
import pysam
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from builtins import range
from future.utils import lrange
from optparse import OptionParser
from random import random, randint, randrange
from scipy.stats import beta, laplace
from numpy.random import multivariate_normal, choice
from pybedtools import BedTool

MAX_TEST = 10000 #number of maximal test iterations, can influence runtime

#classes
class HelpfulOptionParser(OptionParser):
	"""An OptionParser that prints full help on errors."""
	def error(self, msg):
		self.print_help(sys.stderr)
		self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


#helper functions
def _callback_list(option, opt, value, parser): #get a list of options and save it
	try:
		tmplist = list(map(lambda x: int(x), value.split(',')))
		setattr(parser.values, option.dest, tmplist)  
	except Exception as e:
		print("got exception: %r, please check input parameters" % (e,))
		exit()

def _callback_matrix(option, opt, value, parser): #get a matrix of options and save it
	try:
		final_option=[]
		for entry in value.split(';'):
			final_option.append(list(map(lambda x: int(x), entry.split(','))))
		setattr(parser.values, option.dest, final_option)
	except Exception as e:
		print("got exception: %r, please check input parameters" % (e,))
		exit()

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

def _scale_exp( nfrags, read_frac, loc, scale):
	"""Scale number_frags result based on Exponential distribution"""
	top_beta_scale = sp.stats.expon.pdf(loc,loc=loc, scale=scale)#to get to 100% at size loc
	scale_factor=1/top_beta_scale #this is the factor to set beta_frac to 100% at loc
	beta_scale = sp.stats.expon.pdf(nfrags,loc=loc, scale=scale) *scale_factor
	read_frac_old = read_frac
	
	if(beta_scale == 0):
		read_frac = read_frac 
	else:
		read_frac = ((read_frac-0.5)*beta_scale)+0.5 #center and scale it based on beta_scale from exponential distribution, we do this to reduce the read_frac towards 0.5 for high values
	
	#print("nfrags: %s read_frac_in: %s scaling_factor: %s read_frac_out: %s" % (nfrags ,read_frac_old, beta_scale, read_frac))
	return read_frac


#main
if __name__ == '__main__':
	#handle input arguments
	parser = HelpfulOptionParser(usage=__doc__)
	
	parser.add_option("-c", "--chrom", default='chr19', dest="chromosome", type="string", help="Chromosome used for simulation [default: %default]")

	parser.add_option("--read-count-scaling", default="none", dest="frag_count_scaling", type="string", help="Scaling of read distribution, no scaling, scaling of beta result based on read counts (with exp) or scaling of read counts based on beta result (with laplace) : none , frag , beta [default: %default]")
	parser.add_option("--read-count-lp-scale", default=0.1, dest="frag_count_lp_sc", type="float", help="Scale for Laplace distribution if read-count-scaling is frag [default: %default]")
	parser.add_option("--read-count-ex-loc", default=10, dest="frag_count_ex_lo", type="float", help="Loc for exponential distribution if read-count-scaling is beta [default: %default]")
	parser.add_option("--read-count-ex-scale", default=100, dest="frag_count_ex_sc", type="float", help="Scale for exponential distribution if read-count-scaling is beta [default: %default]")
	
	parser.add_option("--beta", default=[0.5, 0.5], dest="beta_values", type="string", action='callback', callback=_callback_list_float, help="Alpha and Beta of Beta-distribution [default: %default]")
	
	(options, args) = parser.parse_args()
	num_args = 4
	
	#read in required input arguments
	if len(args) != num_args:
		parser.error("At minimum %s parameters are requred to run this tool" % num_args)
	input_bam_file = args[0]
	inbed_filename = args[1]
	in_chromsizes = args[2]
	num_replicates = int(args[3])
	
	#get sample bam file...
	input_bam = pysam.AlignmentFile(input_bam_file,'rb')
	subsample_chrom = options.chromosome 
	
	#read bed
	input_bed = BedTool(inbed_filename)
	
	region_count=1
	fragment_counts=[]
	sample_counts=[]
	sumcounts=[]
	
	#for each bed entry get reads...
	sample1_replicate_counts = [0] * num_replicates
	sample2_replicate_counts = [0] * num_replicates
	sample1_replicate_list = [[] for i in lrange(num_replicates)]
	sample2_replicate_list = [[] for i in lrange(num_replicates)]
	multi_factor = 2*num_replicates #per read multiply its use, so that we can definitely fill 2 samples with x replicates
	
	for region in input_bed:
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
		elif (options.frag_count_scaling == "frag") :#nfrags scaling via lognorm distribution based on fragment counts, number_frags not changed
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
				else:
					sample2_count += 1
			else:
				continue
			break
			
		#save result for this protein
		sample_counts.append([sample1_count,sample2_count])
		sumcounts.append(number_frags)

	#MA-plot
	xma=[]
	yma=[]
	counts1=[]
	counts2=[]
	for entry in sample_counts:
		counts1.append(entry[0])
		counts2.append(entry[1])
		if(entry[1] > 0 and entry[0] > 0):
			x=math.log(entry[0]+entry[1]/2,2)#log2 mean
			xma.append(x) 
			y=math.log(entry[0]/entry[1],2)#log2 fold-change
			yma.append(y) 
	
	
	print("max: sum %s s1 %s s2 %s" % (max(sumcounts),max(counts1), max(counts2)))
	
	#hist and ma-plot
	fixed_range=[0,200]
	fixed_bins=500
	fig, axs = plt.subplots(4)
	
	axs[0].set_title('fragments all samples')
	axs[0].hist(sumcounts,density=False, bins=fixed_bins, range=fixed_range)
	axs[1].set_title('fragments sample1')
	axs[1].hist(counts1,density=False, bins=fixed_bins, range=fixed_range)
	axs[2].set_title('fragments sample2')
	axs[2].hist(counts2,density=False, bins=fixed_bins, range=fixed_range)
	axs[3].set_title('MA-plot')
	axs[3].scatter(xma, yma)
	
	plt.subplots_adjust(hspace = 0.4)
	plt.show()
	
	#plt.scatter(xma, yma)
	#plt.show()
