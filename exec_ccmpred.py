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
import re

nthreads = 30
ccmpred_exe = 'ccmpred'

"""
Run CCMpred and get the normalized contact matrix
"""
def exec_ccmpred(inputdir, outputdir, alnfile):
	ccmpred_path = os.environ['CCMPRED']
	outfile = outputdir + alnfile + '.ccmpred.out'
	if not ccmpred_path.endswith('/'):
		ccmpred_path += '/' 
	ccmpred_cmd = ccmpred_path + ccmpred_exe + ' -t ' + str(nthreads) + ' -R ' + inputdir + alnfile + ' ' + outfile 
	errstr = ''
	rv = os.system(ccmpred_cmd + ' >/dev/null 2>&1')
	if rv:
		errstr = 'ERROR encountered while running CCMpred'
		return rv, errstr
	else:
		return rv, outfile
