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



""" 
The master script for KFC-E.
Input arguments:
	inputdir: Input directory having the PDB complex_file and alignment file [Required]
	complex_file: Name of the PDB file [Required]. Must contain only two chains, the ones which will be docked [Required]
	alignment_file: Name of the file having the paired alignments [Required]
	outputdir: The name of the output directory where all the results will be placed [Required]
	kfc_e_outfile: The name of the file (in the output directory) into which the KFC-E output will be written to [Optional]
				   If no file name is provided, the script will write the output into the default file				

This script will perform the following tasks
	1. Process the PDB file by renumbering the residues serially, remove hydrogens and any ligands, 
	   include only the 'A' locations (in case of alternative locations) and write the processed coordinates 
	   into a separate file. Rename the chains to A and B.
	2. Run CCMPred calculations on the alignment file
	3. Run ZDOCK to exhaustively sample the binding modes
	4. Run KFC2 predictions on each sampled binding pose from ZDOCK
	5. Run calculations for interface contacts (at 7 Ang cutoff) for each sampled pose
	6. Run calculations for residue conservation
	6. Score each pose with KFC-E

Required environment variables:
	ZDOCK: Set to the path of the ZDOCK executable
	RATE4SITE; Set to the path of the RATE4SITE executable
	CCMPRED: Set to the path of the ccmpred executable
	KFC2: Set to the path of the KFC2 executable
	KFCE_HOME: Set to the directory with KFC_E codebase
	CLUSTALO: Set to the path containing the clustal omega executable
"""

import os
import sys
import argparse
import pdbutils
import shutil
import exec_ccmpred as exec_ccmpred 
import exec_zdock as exec_zdock
import exec_kfc2 as exec_kfc2
import calc_contacts as calc_contacts
import calc_conservation as cons
import KFC_E_score as KFC_E_score

"""
Checks whether all the required environment variables are set 
or not.

Input Arguments: None

Returns: 
	a. Error string if all variables are not set
	b. Null, if all variables are set
"""
def check_env():
	req_env_vars = ['ZDOCK', 'RATE4SITE', 'CCMPRED', 'KFC2', 'KFCE_HOME', 'CLUSTALO']
	unset_vars = []
	for var_i in req_env_vars:
		if not var_i in os.environ.keys():
			unset_vars.append(var_i)
	return unset_vars

"""
Checks the input PDB and alignment files
Input Arguments:
	inputdir - The input directory
	pdbfile - The pdb file 
	alnfile - The paired alignment file

Returns: 
	a. Error string if files are not OK
"""
def checkinput(inputdir, pdbfile, alnfile):
	outstr = ''
	if not os.path.isdir(inputdir):
		outstr += 'No such directory ' + inputdir + '\n'
	if not os.path.isfile(inputdir+pdbfile):
		outstr += 'No such file ' + pdbfile + '\n'
	if not os.path.isfile(inputdir+alnfile):
		outstr += 'No such file ' + alnfile + '\n'
	return outstr

"""
Checks the total number of chains in the pdb file
Input Arguments:
	inputdir - The input directory
	pdbfile - The name of the pdb file

Returns: The number of chains

"""
def get_numchains(inputdir, pdbfile):
	pstruct = pdbutils.parse_pdb_as_dict(inputdir,pdbfile,[],[])	
	return len(pstruct['CHAINORDER'])

"""
Makes sure that the PDB file is parsed properly into a dictionary 
having the same number of [X,Y,Z] coordinates, residue names,
residue ids, atom names and atom numbers.

Input: The pdb file parsed into a dictionary

Returns: 
	rv - A boolean variable that tells whether the dictionary meets the minimum 
		requirements
	errstr - The error string with the error message
"""
def checkpdbdict(pdbstruct):
	chainorder = pdbstruct['CHAINORDER']
	rv = 0
	errstr = ''
	for chain in chainorder:
		check_keys = ['COORD', 'RESNAME', 'RESID', 'ATOMNAME', 'ATOMNUM']
		numelemlist = []
		for key_i in check_keys:
			numelemlist.append(len(pdbstruct[chain][key_i]))
		if len(set(numelemlist)) > 1:
			rv = 1
			errstr += 'Different number records for chain ' + str(chain) + '\n'
	return rv, errstr	

"""
The main caller function
"""
def main(args):
	inputdir = args['inputdir']
	pdbfile = args['pdbfile']
	alnfile = args['alnfile']
	outputdir = args['outputdir']
	
	# Preliminary checks
	if args['nprocs'] == None: # If number of processes not assigned, set it to 1
		nprocs = 1
	else:
		nprocs = args['nprocs']
	if not inputdir.endswith('/'):
		inputdir += '/'
	if not outputdir.endswith('/'):
		outputdir += '/'
	if not pdbfile.endswith('.pdb'):
		sys.exit('\nERROR! PDB file must end with .pdb')
		
	# Check whether the input files are present or not
	print("Checking input files...", end='', flush=True)
	rv = checkinput(inputdir, pdbfile, alnfile)
	if rv:
		sys.exit('\nERROR! ' + rv)
	print("OK!")

	# Perform a check for the required environment variables
	print("Checking environment variables...", end='', flush=True)
	rv = check_env()
	if len(rv) > 0:
		rv = ", ".join(rv)
		sys.exit("\nERROR! Required environment variable(s) not set: " + rv)
	print("OK!")

	# Check if PDB file has more than two chains
	print("Checking number of chains in PDB file...", end="", flush=True)
	numch = get_numchains(inputdir, pdbfile)
	if numch > 2:
		sys.exit("\nERROR! PDB file has more than 2 chains! Must have only 2 chains!")
	print("OK!")

	# Create the output directory
	if not os.path.isdir(outputdir):
		print("Creating output directory...", end="", flush=True)
		try:
			os.mkdir(outputdir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
		print("Done!")
	
	# Get the input and output directories with their absolute paths
	workdir = os.getcwd()
	os.chdir(inputdir)
	inputdir = os.getcwd()+'/'
	os.chdir(workdir)
	os.chdir(outputdir)
	outputdir = os.getcwd()+'/'
	os.chdir(workdir)

	# Parse the PDB file as a dictionary
	skiph = 1 # Skip hydrogen atoms in the PDB
	print("Parsing and processing input pdb file...", end="", flush=True)
	newchains = ['A', 'B']
	pdbstruct = pdbutils.parse_pdb_as_dict(inputdir, pdbfile, skiph, newchains)
	# Check to make sure that all the fields of pdbstruct have the same number
	# of elements
	rv, errstr = checkpdbdict(pdbstruct)
	if rv == 1:
		sys.exit("\nERROR parsing PDB file into dictionary: " + errstr)
	# Re-number the PDB records with atom numbers and residue numbers starting with 1
	pdbstruct_new = pdbutils.renumberpdb(pdbstruct)
	pdbutils.writepdb(pdbstruct_new, outputdir, pdbfile, [])
	print("Done!")


	# Write the individual chains into different PDB files
	chainA_pdb = pdbfile.replace('.pdb', '_A.pdb')
	chainB_pdb = pdbfile.replace('.pdb', '_B.pdb')
	pdbutils.writepdb(pdbstruct_new, outputdir, chainA_pdb, 'A')
	pdbutils.writepdb(pdbstruct_new, outputdir, chainB_pdb, 'B')

	# Run CCMPred calculations
	ccmpred_outputfile = outputdir + alnfile + '.ccmpred.out'
#	print ("Running CCMpred calculations...", end='', flush=True)
#	rv, msg = exec_ccmpred.exec_ccmpred(inputdir, outputdir, alnfile)
#	if rv:
#		sys.exit('\nERROR in CCMpred execution')

	min_files = 54000
	
	# Run ZDOCK
	pdbprefix = pdbfile.replace('.pdb', '')
	zdockoutputdir = outputdir + pdbprefix + '_zdock_out/'
	if not os.path.isdir(zdockoutputdir):
		try:
			os.mkdir(zdockoutputdir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	print ("Running ZDOCK...", end='', flush=True)
	rv, numfiles = exec_zdock.exec_zdock(chainA_pdb,chainB_pdb,outputdir,zdockoutputdir)
	if rv and (numfiles < min_files):
		sys.exit('\nERROR running ZDOCK')
	elif not rv and (numfiles < min_files):
		sys.exit('\nERROR with ZDOCK! Not all poses generated')
	else:
		print ("\nWrote ZDOCK sampled poses to ",zdockoutputdir)

	# Perform hotspot calculations with KFC2
	kfc2outputdir = outputdir + pdbprefix + '_kfc2_hotspots/'
	if os.path.isdir(kfc2outputdir):
		shutil.rmtree(kfc2outputdir)
	try:
		os.mkdir(kfc2outputdir)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass
	print ("Running KFC2 hotspot calculations for docked poses...", end='', flush=True)
	numfiles = exec_kfc2.exec_kfc2(zdockoutputdir, kfc2outputdir, nprocs)
	if numfiles < min_files:
		sys.exit('\nERROR with KFC2 calculations. Minimum number of KFC2 files not found')
	else:
		print ("Done!")
	
	# Perform contact calculations for docked poses
	contacts_outputdir = outputdir + pdbprefix + '_contacts/'
	if os.path.isdir(contacts_outputdir):
		shutil.rmtree(contacts_outputdir)
	try:
		os.mkdir(contacts_outputdir)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass
	print ("Calculating contacts...", end='', flush=True)
	cutoff_l = 4
	cutoff_u = 7
	numfiles = calc_contacts.calc_contacts(zdockoutputdir, contacts_outputdir, cutoff_l, cutoff_u, nprocs)
	if numfiles < min_files:
		sys.exit("\nERROR with contact calculations. Required number of contact files (", min_files, ") not found")
	else:
		print ("Done!")

	# Perform calculations for residue conservation
	print ("Running calculations for residue conservation...", end='', flush=True)
	cons_dir = outputdir + pdbprefix + '_cons_scores/'
	if os.path.isdir(cons_dir):
		shutil.rmtree(cons_dir)
	try:
		os.mkdir(cons_dir)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass	
	N,cumul_cons_file = cons.conservation(alnfile, inputdir, cons_dir)
	
	# Perform scoring with KFC-E
	scorefile, scorematrix = KFC_E_score.scoreposes_wrapper(outputdir+pdbfile, inputdir+alnfile, ccmpred_outputfile, cumul_cons_file, zdockoutputdir, contacts_outputdir, kfc2outputdir, outputdir)
	if len(scorematrix) < min_files:
		sys.exit('\nERROR in scoring poses!')
	else:
		print ("Done!\nWrote KFC-E scores to ", scorefile)


"""
Caller to the main function
"""
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--inputdir", type=str, required=True, help='Input directory containing the PDB and the alignment files [REQUIRED]')
	parser.add_argument("--pdbfile", type=str, required=True, help='Protein complex PDB file [REQUIRED]')
	parser.add_argument("--alnfile", type=str, required=True, help='Alignment file in ccmpred format [REQUIRED]')
	parser.add_argument("--outputdir", type=str, required=True, help='Output directory [REQUIRED]')
	parser.add_argument("--nprocs", type=int, required=False, help='Number of processes [OPTIONAL]')

	args = vars(parser.parse_args())
	main(args)
