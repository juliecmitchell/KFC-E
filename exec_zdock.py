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
import pdbutils as ptils
from shutil import copyfile
import glob

"""
Decide on which chain is the receptor and which the ligand
We will consider the larger chain (the one with more residues)
as the receptor and the smaller chain as the ligand
"""
def assign_rec_lig_chains(outputdir, chain_1, chain_2):
	chain_1 = outputdir + chain_1
	chain_2 = outputdir + chain_2
	parsed_chain1 = ptils.parse_coordinates(chain_1)
	parsed_chain2 = ptils.parse_coordinates(chain_2)
	nres1 = len(parsed_chain1['UNIQ_RES_ID'])
	nres2 = len(parsed_chain2['UNIQ_RES_ID'])
	if nres1 > nres2:
		return ('A','B')
	else:
		return ('B','A')

"""
Run Z dock for a given pdb file
"""
def run_z_dock(pdb_prefix, rec_chain, lig_chain, outputdir, zdockoutputdir):
	zdock_path = os.environ['ZDOCK']
	N_poses = 80000
	#N_poses = 10
	dense_sampling_flag = 1
	
	# Copy the ligand and receptor files to the zdockoutput dir
	rec_file = pdb_prefix + '_' + rec_chain + '.pdb'
	lig_file = pdb_prefix + '_' + lig_chain + '.pdb'
	copyfile(outputdir+rec_file, zdockoutputdir+rec_file)
	copyfile(outputdir+lig_file, zdockoutputdir+lig_file)

	# Copy the required files for running Z dock
	if not zdock_path.endswith('/'):
		zdock_path += '/'
	copyfile(zdock_path+'create_lig', zdockoutputdir+'create_lig')
	copyfile(zdock_path+'uniCHARMM', zdockoutputdir+'uniCHARMM')
	os.chmod(zdockoutputdir+'create_lig', 751) # provide the permissions for user (rwx), others (execute), group (read and execute)	
	
	# Get the current working directory
	work_dir = os.getcwd()
	# Change directory to the output dir
	os.chdir(zdockoutputdir)

	# Mark the surface of the receptor and ligand
	rec_m_file = pdb_prefix + '_' + rec_chain + '_m.pdb'
	lig_m_file = pdb_prefix + '_' + lig_chain + '_m.pdb'
	cmd1 = zdock_path + 'mark_sur' + ' ' + rec_file + ' ' + rec_m_file
	cmd2 = zdock_path + 'mark_sur' + ' ' + lig_file + ' ' + lig_m_file
	#print ("Marking surface for receptor... ", flush=True, end='')
	os.system(cmd1)
	#print ("Done!")
	#print ("Marking surface for ligand... ",flush=True, end='')
	os.system(cmd2)
	#print ("Done!")

	# Run Z-dock
	cmd3 = zdock_path + 'zdock'
	if dense_sampling_flag == 1:
		cmd3 = cmd3 +  ' -D'
	
	z_dock_output_file = pdb_prefix + '_AB.zdock.out'
	cmd3 = cmd3 + ' -N ' + str(N_poses) + ' -R ' +  rec_m_file + ' -L ' + lig_m_file + ' -o ' + z_dock_output_file
	#print ("cmd3 = ",cmd3)
	#print ("Running Z-dock... ", end='', flush=True)
	os.system(cmd3)
	#print ("Done!")

	# Translate Z-dock output into PDB poses
	cmd4 = 'perl ' + zdock_path + 'create_2.pl' +  ' ' + z_dock_output_file + ' ' + str(N_poses) + ' ' + pdb_prefix + '_AB'
	#print ("Getting docked poses... ", end='', flush=True)
	#print ("cmd4 = ", cmd4)
	rv = os.system(cmd4)

	# Clean up by removing the uniCHARMM and create_lig files that were copied and also the receptor_m, ligand_m files, receptor and ligand PDB files.
	#print ("Cleaning up ...", end='', flush=True)
	os.remove('create_lig')
	os.remove('uniCHARMM')
	os.remove(rec_m_file)
	os.remove(lig_m_file)
	if os.path.isfile(rec_file):
		os.remove(rec_file)
	if os.path.isfile(lig_file):
		os.remove(lig_file)
	#print("Done!\n")

	numfiles = len(glob.glob('*.pdb')) # Check for the total number of files. Should be 54,000 when using exhaustive sampling in ZDOCK.
	os.chdir(work_dir) # Return to working directory
	return rv, numfiles

"""
Caller function for docking
"""
def exec_zdock(chain1_pdb, chain2_pdb, outputdir, zdockoutputdir):
	rec_chain, lig_chain = assign_rec_lig_chains(outputdir, chain1_pdb, chain2_pdb)
	pdb_prefix = chain1_pdb.replace('_A.pdb', '')
	return run_z_dock(pdb_prefix, rec_chain, lig_chain, outputdir, zdockoutputdir)
		
