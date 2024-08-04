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
import re
import sys 
import pdbutils as ptils
from shutil import copyfile
from multiprocessing import Pool
import glob
import shutil

"""
Run the kfc2 predictions for each file 
"""
def run_kfc2(listiterator_i):
	inputdir = listiterator_i[0]	
	outputdir = listiterator_i[1]
	pdbfile = listiterator_i[2]
	pdb_base_name = pdbfile.replace('.pdb', '')
	chain_1 = 'A'
	chain_2 = 'B'
	
	# Get the full path for the input and output directories (in case the user has not provided the full path)
	workdir = os.getcwd()
	os.chdir(inputdir)
	inputdir_fullpath = os.getcwd()
	if not inputdir_fullpath.endswith('/'):
		inputdir_fullpath += '/'
	os.chdir(workdir)
	os.chdir(outputdir)
	outputdir_fullpath = os.getcwd()
	if not outputdir_fullpath.endswith('/'):
		outputdir_fullpath += '/'
	
	kfc_output_file = pdb_base_name+'.kfc.results' # This is the file we are interested in! It will contain both the KFC A and KFC B predictions.
	if os.path.isfile(outputdir_fullpath + kfc_output_file):
		print ("\nFile ", kfc_output_file, " already present. Skipping calculations for this pose!")
		return
	# Create a temporary directory inside outptdir2. Here we will store all the output for kfc2 and then, retain only the
	# hot spot prediction information for this pose.
	tmpdir = pdb_base_name +  '/'
	if not os.path.isdir(tmpdir):
		try:
			os.mkdir(tmpdir)
		except:
			print ("Cannot create directory!\n")
			return
	try:
		copyfile(inputdir_fullpath + pdbfile, outputdir_fullpath + tmpdir + pdbfile) # copy the pose_i to the new directory
	except IOError as e:
		print ("Unable to copy file.", e)
		return
	os.chdir(tmpdir)

	# Let's get the predictions from KFC2 for the given PDB file
	kfc_path = os.environ['KFC2']
	if not kfc_path.endswith('/'):
		kfc_path += '/'
	cmd = kfc_path + 'kfc2.sh' + ' ' + pdbfile + ' ' + chain_1 + ' ' + chain_2 + ' ' + pdb_base_name
	#print ("cmd = ", cmd)
	os.system(cmd + '>/dev/null 2>&1')
	#os.system(cmd)
	#kfc_output_file = pdb_base_name+'.kfc.results' # This is the file we are interested in! It will contain both the KFC A and KFC B predictions.

	copyfile(outputdir_fullpath+tmpdir+kfc_output_file, outputdir_fullpath + kfc_output_file)
	os.chdir(workdir)
	shutil.rmtree(outputdir_fullpath+tmpdir, ignore_errors=True) # remove the temporary directory

"""
Caller function for run_kfc2
"""
def exec_kfc2(inputdir, outputdir, nprocs):
	workdir = os.getcwd()
	os.chdir(inputdir)
	pdbfiles = glob.glob('*.pdb')
	os.chdir(workdir)
	list_iterator = [[inputdir,outputdir,file_i] for file_i in pdbfiles]
	try:
		p = Pool(nprocs)
		p.map(run_kfc2, list_iterator)
	finally:
		p.close()
		p.join()
	os.chdir(outputdir)
	numfiles = len(glob.glob('*.kfc.results'))
	os.chdir(workdir)
	return numfiles
