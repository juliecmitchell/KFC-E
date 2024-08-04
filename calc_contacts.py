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
from multiprocessing import Pool
import pdbutils as ptils
import glob
import numpy as np

"""
Calculate the contacts in the interface of the given pdb file
"""
def contacts(list_iter):
	pdb_i, inputdir, outputdir, cutoff_l, cutoff_u = list_iter
	pdb_i_basename = pdb_i.replace('.pdb', '')
	pdb_i_fq = inputdir + pdb_i
	outfile1 = outputdir + pdb_i_basename + '_A.pdb'
	outfile2 = outputdir + pdb_i_basename + '_B.pdb'
	ptils.split_pdb_2_chains(pdb_i_fq, 'A', 'B', outfile1, outfile2)
	os.chdir(outputdir)
	file_1 = pdb_i_basename + '_A.pdb'
	file_2 = pdb_i_basename + '_B.pdb'
	kfce_home = os.environ['KFCE_HOME']
	if not kfce_home.endswith('/'):
		kfce_home += '/'
	contact_prg = kfce_home + 'contacts'
	cmd = contact_prg + ' ' + file_1 + ' ' + file_2 + ' ' + str(cutoff_l) + ' ' + str(cutoff_u)
	os.system(cmd)
	
	# Remove the split pdb files and also the coresponding *.fuc files. The .fuc files
	# sometimes do not give accurate contacts. We will only retain the .fc file and then 
	# calculate the unique contacts from this file
	fuc_file = file_2 + '.fuc'
	os.remove(file_1)
	os.remove(file_2)
	os.remove(fuc_file)	
	fc_file = file_2 + '.fc'
	new_fuc_file = pdb_i_basename + '.fuc'
	new_fc_file = pdb_i_basename + '.fc'
	os.rename(fc_file, new_fc_file) # rename the fc file
	
	# Calculate the unique contacts
	uniq_cont = []
	fh = open(new_fc_file, 'r')
	for line_i in fh: 
		r1,r2 = int(line_i.split()[1]), int(line_i.split()[6])
		if [r1,r2] not in uniq_cont and [r2,r1] not in uniq_cont:
			uniq_cont.append([r1,r2])
	fh.close()	
	np.savetxt(new_fuc_file, uniq_cont, delimiter=',', fmt='%d')

"""
Caller for the contacts function
"""
def calc_contacts(inputdir, outputdir, cutoff_l, cutoff_u, nproc):
	workdir = os.getcwd()
	os.chdir(inputdir)
	all_pdbs = glob.glob('*.pdb')
	os.chdir(workdir)
	list_iter = [[pdb_i, inputdir, outputdir, cutoff_l, cutoff_u] for pdb_i in all_pdbs]
	try:
		p = Pool(nproc)
		p.map(contacts, list_iter)
	finally:
		p.close()
		p.join()
	os.chdir(outputdir)
	numfiles = len(glob.glob('*.fc'))
	os.chdir(inputdir)
	return numfiles

