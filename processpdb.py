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
import numpy
import sys

"""
This function will parse a given PDB file and retain only the ATOM records and eliminate all the HETATM records
Also, in case of alternate locations like ASER and BSER, it will retain the 'A' records only.
"""
def process_pdb(pdbfile,inputdir):
	for pdbfile in pdbfiles:
		if not re.match('.*?\.pdb$', pdbfile):
			continue
		else:
			pdbfile_fq = inputdir  + pdbfile
			outfile_fq = outputdir + pdbfile
			fh1 = open(pdbfile_fq, 'r')
			fh2 = open(outfile_fq, 'w')
			print("Processing file ", pdbfile)
			for line in fh1:
				if re.match('^ATOM', line):
					alt_loc = line[16] # alternate location if any
					alt_loc = alt_loc.replace(' ','')
					if alt_loc == 'A': # we will consider by default the 'A' location
						fh2.write(line)
					elif alt_loc == '':
						fh2.write(line)
					else:
						continue
						
				elif re.match('^TER', line):
					fh2.write(line)
				else:
					continue
			fh1.close()
			fh2.close()
			print("Done!\n")

