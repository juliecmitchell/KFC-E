# MIT License
#
# Copyright (c) 2019 Julie C. Mitchell and Oak Ridge National Laboratory
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


import os
import sys
import numpy as np
from multiprocessing import Pool
import random
from shutil import copyfile
import glob
import re

"""
Create an input alignment file in the temporary alignment directory for rate4site
"""
def create_alignment(aln_subset_file, sampled_aln):
	fh = open(aln_subset_file, 'w')
	counter = 0
	for s in sampled_aln:
		fh.write(">" + str(counter) + "\n" + s + "\n")
		counter += 1
	fh.close()

"""
Decide on the sample size and number of samples of 
sequences to be considered from the alignment file
"""
def get_aln_subset_size(alnfile):
	numseqs = len(np.loadtxt(alnfile, dtype=str))
	if numseqs <= 120:
		N_iter = 1
		N_max = numseqs
	else:
		N_max = 120
		N_iter = int(np.ceil(numseqs/N_max)*5)
	return N_max, N_iter	

"""
Calculate the average conservation score from the sampled alignments
"""
def calc_cumulative_cons(cumul_cons_file):
	cons_files = glob.glob('*.cons')
	all_scores = []
	if os.path.isfile(cumul_cons_file): # Nothing to do if the cumulative conservation file is already present
		print ("Conservation file already present. Nothing to do!")
		return
	for cons_file_i in cons_files:
		scores_i = []
		fh = open(cons_file_i, 'r')
		for line in fh:
			if line.startswith('#') or re.match('^\s+?$', line): # skip lines starting with # or blank lines
				continue
			else:
				line_split = line.split() # split by white space and then append score
				scores_i.append(float(line_split[2]))
		fh.close()
		all_scores.append(scores_i)
	# Calculate the median scores (we will stick to median since it is a more robust estimator)
	median_scores = np.median(np.array(all_scores), axis=0)
	np.savetxt(cumul_cons_file, median_scores, fmt='%0.4f')

"""
Wrapper function for calculating the conservation scores
"""
def calc_cons(alnfile, inputdir, outputdir):
	copyfile(inputdir + alnfile, outputdir + alnfile)
	workdir = os.getcwd()
	os.chdir(outputdir)
	N_max,N_iter = get_aln_subset_size(alnfile)
	aln_seqs = [line.rstrip() for line in open(alnfile, 'r')]
	alnfile_prefix = alnfile.rsplit('.', 1)[0]
	rate4site_exec = os.environ['RATE4SITE']
	if not rate4site_exec.endswith('/'):
		rate4site_exec += '/'
	rate4site_exec += 'rate4site'
	# Calculate conservation for N_max subset of randomly selected sequences N_iter
	# times.
	for i in range(0,N_iter):
		aln_subset_file = alnfile_prefix + '.' + str(i) + '.aln'
		cons_subset_file = alnfile_prefix + '.' + str(i) + '.cons'
		try:
			# Sample N_max random sequences
			sampled_ind = random.sample(list(range(1,len(aln_seqs))), N_max)
			sampled_ind = [0] + sampled_ind # include the 0th index which is the native sequence
			sampled_aln = [aln_seqs[i] for i in sampled_ind]

			# Create the alignment file for the sampled subset of sequences
			create_alignment(aln_subset_file, sampled_aln)

			# Run rate4site on the alignment file
			cmd = rate4site_exec + ' -s ' + aln_subset_file +  ' -o ' + cons_subset_file + ' >/dev/null 2>&1'
			#print ("Running rate4site for ", alnfile_prefix, " , iter ", i)
			os.system(cmd)
			#print ("Done!")
		except Exception as e: # Do not break the process in case of an exception
			print ("Encountered exception in rate4site calculations for ", aln_subset_file, ". Will skip!!! Exception : ", e)
			continue
	cumul_cons_file = outputdir + alnfile_prefix + '.' + 'cumulative.cons'
	calc_cumulative_cons(cumul_cons_file)
	os.chdir(workdir)
	return cumul_cons_file

"""
Caller to the wrapper function for calculating conservation
"""
def conservation(alnfile, inputdir, outputdir):
	cumul_cons_file = calc_cons(alnfile, inputdir, outputdir)
	N = len(np.loadtxt(cumul_cons_file, dtype=float))
	return N, cumul_cons_file
	
