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


# Set of functions for parsing PDB files

import numpy as np
import re
import sys
import copy

"""
Parse the coordinates from the PDB file
"""
def parse_coordinates(filename):
	# parse coordinate information from a PDB file and return information
	# as a dictionary
	coords = [] # multi dimensional list to hold the coordinate information for each residue
	all_residues = [] # store the residue names
	all_res_ids = [] # store the residue ids
	all_res_chain_ids = []
	unique_residues = [] # store the unique residue names
	unique_residue_ids = [] # store the unique residue ids
	unique_res_chain_ids = [] # store the chain-residue identifier (e.g., A12,A13, B13, B15 and so on where A and B are the chain IDs)
	parsed_info = {} # dictionary to store the parsed information
	fh = open (filename, 'r')
	for line in fh:
		if re.match('^ATOM',line):
			line = line.rstrip()
			chain = line[21]
			X_coord = line[30:38]
			Y_coord = line[38:46]
			Z_coord = line[46:54]
			
			# Strip white space
			X_coord = re.sub('\s+','',X_coord)
			Y_coord = re.sub('\s+','',Y_coord)
			Z_coord = re.sub('\s+','',Z_coord)
			
			# convert from string to numeric
			X_coord = float(X_coord)
			Y_coord = float(Y_coord) 
			Z_coord = float(Z_coord) 
			coords.append([X_coord, Y_coord, Z_coord])						
			res_name = line[17:20]
			res_id = line[22:26]
			
			res_name = re.sub('\s+','',res_name)
			res_id = re.sub('\s+','',res_id)
			res_chain_id = chain+res_id

			all_res_chain_ids.append(res_chain_id)
			res_id = int(res_id)
			
			# include all residue names and ids
			all_residues.append(res_name) 
			all_res_ids.append(res_id)
			
			# check for unique residue ids i.e., remove redundancy occuring from side chain
			if len(unique_residues) == 0: # if empty then include current residue name and id
				unique_residues.append(res_name)
				unique_residue_ids.append(res_id)
				unique_res_chain_ids.append(res_chain_id)
			else: # Include the residue id, name and chain-residue identfier if previous id in the list 
					# is different from the current one. 
				prev_uniq_res_id = unique_residue_ids[-1]
				if res_id != prev_uniq_res_id:
					unique_residues.append(res_name)
					unique_residue_ids.append(res_id)
					unique_res_chain_ids.append(res_chain_id)
		else:
			#pass
			continue
	coord_np = np.array(coords)
	
	# Store the parsed information in a dictionary and return to caller
	parsed_info["COORD"] = coord_np
	parsed_info["UNIQ_RES"] = np.array(unique_residues)
	parsed_info["UNIQ_RES_ID"] = np.array(unique_residue_ids)
	parsed_info["ALL_RES"] = np.array(all_residues)
	parsed_info["ALL_RES_ID"] = np.array(all_res_ids)
	parsed_info["UNIQ_RES_CHAIN_ID"] = np.array(unique_res_chain_ids)
	parsed_info["ALL_RES_CHAIN_ID"] = np.array(all_res_chain_ids)
	return parsed_info

"""
Parse the PDB file as a dictionary

Input arguments:
	inputdir - Name of the input directory
	filename - Name of the PDB file in the input directory to be parse
	skiph - Boolean variable. When set to 1, skips H atoms otherwise not.
	newchains - List of new chains into which we will rename the existing chains
"""
# parse the PDB file as a dictionary 
def parse_pdb_as_dict(inputdir, filename, skiph, newchains):
	# parse the contents of a PDB file as a dictionary. This is especially useful
	# if the PDB file has more than one chain.
	PDB_struct = {}
	filename = inputdir + filename
	fh = open(filename, 'r')
	chainorder = [] # The original ordering of chains in the given PDB file
	chainind = -1
	#print(newchains)
	for line in fh:
		if re.match('^ATOM', line):
			chain = line[21]
			#print (chain)
			# The chain to include in the parsed dictionary
			# will depend on whether newchains is set or not.

			if chain not in chainorder:
				chainind = chainind + 1
				#print(chainind)
				chainorder.append(chain)
			if newchains:
				chain = newchains[chainind]

			X_coord = line[30:38]
			Y_coord = line[38:46]
			Z_coord = line[46:54]
			
			# Strip white space
			X_coord = re.sub('\s+','',X_coord)
			Y_coord = re.sub('\s+','',Y_coord)
			Z_coord = re.sub('\s+','',Z_coord)
			
			# convert from string to numeric
			X_coord = float(X_coord)
			Y_coord = float(Y_coord) 
			Z_coord = float(Z_coord)	
			
			res_name = line[17:20]
			res_id = line[22:26]
			
			res_name = re.sub('\s+','',res_name)
			res_id = re.sub('\s+','',res_id)
			res_id = int(res_id)
			atomnum = int(line[6:11].replace(' ', ''))
			atomname = line[13:16].replace(' ','')
			if skiph and (atomname.startswith('H') or (atomname[0].isdigit() and atomname[1] == 'H')): # skip hydrogen information if skiph is set
				continue
			altloc = line[16].replace(' ','')
			if altloc and altloc != 'A': # skip if the alternate location is set and is not A
				continue
			occ = float(line[54:60].replace(' ', '')) # occupancy
			bfactor = float(line[60:66].replace(' ', '')) # bfactor
			elesym = line[76:78] 
			if chain not in PDB_struct.keys():
				PDB_struct[chain] = dict()
				PDB_struct[chain]['COORD'] = []
				PDB_struct[chain]['RESNAME'] = []
				PDB_struct[chain]['RESID'] = []
				PDB_struct[chain]['ATOMNAME'] = []
				PDB_struct[chain]['ATOMNUM'] = []
				PDB_struct[chain]['ALTLOC'] = []
				PDB_struct[chain]['OCC'] = []
				PDB_struct[chain]['BFACTOR'] = []
				PDB_struct[chain]['ELESYM'] = []
			PDB_struct[chain]['COORD'].append([X_coord, Y_coord, Z_coord])
			PDB_struct[chain]['RESNAME'].append(res_name)
			PDB_struct[chain]['RESID'].append(res_id)
			PDB_struct[chain]['ATOMNUM'].append(atomnum)
			PDB_struct[chain]['ATOMNAME'].append(atomname)
			PDB_struct[chain]['ALTLOC'].append(altloc)
			PDB_struct[chain]['OCC'].append(occ)
			PDB_struct[chain]['BFACTOR'].append(bfactor)
			PDB_struct[chain]['ELESYM'].append(elesym)
	fh.close()
	if newchains:
		chainorder = newchains
	PDB_struct['CHAINORDER'] = chainorder
	return PDB_struct

"""
For a PDB file having 2 chains, write each chain into a separate file.
"""
def split_pdb_2_chains(inputfile, chain1, chain2, outfile1, outfile2):
	fh = open(inputfile, 'r')
	fh1 = open(outfile1, 'w')
	fh2 = open(outfile2, 'w')
	for line in fh:
		if line[21] == chain1:
			fh1.write(line)
		elif line[21] == chain2:
			fh2.write(line)
	fh.close()
	fh1.close()
	fh2.close()

"""
Write the PDB datastructure into a PDB formatted file

"""
def writepdb(pdbstruct, outdir, outfile, chain):
	if not chain:
		chains = pdbstruct['CHAINORDER']
	else:
		chains = chain
	fh = open(outdir+outfile, 'w')
	for chain in chains:
		for i in range(0,len(pdbstruct[chain]['ATOMNAME'])):
			f1 = 'ATOM'
			f2 = pdbstruct[chain]['ATOMNUM'][i]
			f3 = pdbstruct[chain]['ATOMNAME'][i]
			f4 = pdbstruct[chain]['ALTLOC'][i]
			f5 = pdbstruct[chain]['RESNAME'][i]
			f6 = chain
			f7 = pdbstruct[chain]['RESID'][i]
			f8 = '' # code for insertion of residues
			f9 = pdbstruct[chain]['COORD'][i][0] # X coord
			f10 = pdbstruct[chain]['COORD'][i][1] # Y coord
			f11 = pdbstruct[chain]['COORD'][i][2] # Z coord
			f12 = pdbstruct[chain]['OCC'][i]
			f13 = pdbstruct[chain]['BFACTOR'][i]
			f14 = pdbstruct[chain]['ELESYM'][i]
			fmt_str = "{:6s}{:5d}  {:<4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format(f1,f2,f3,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14)
			#fmt_str ="{:6s}{:5d}  {:<4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}"
			#fmt_str = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(f1,f2,f3,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14)
			fh.write(fmt_str + '\n')
	fh.close()


"""
Renumber the residues and atoms in the given PDB structure

Input:
	pdbstruct: Parsed dictionary of the PDB file

Returns:
	pdbstruct_new: New dictionary with updated number of the residue and atoms

"""
def renumberpdb(pdbstruct):
	# Copy the pdb structure into a new dictionary. Simple assignment to a new variable 
	# will not work as modifying the new dictionary will also change the original one. 
	pdbstruct_new = copy.deepcopy(pdbstruct)
	atomnum = 1
	resnum = 1
	prevresnum = resnum
#	print ("pdb struct = ", pdbstruct)
	for chain in pdbstruct['CHAINORDER']:
		pdbstruct_new[chain]['ATOMNUM'] = [] # Reset the atom numbers
		pdbstruct_new[chain]['RESID'] = [] # Reset the res numbers
		for i in range(0,len(pdbstruct[chain]['ATOMNAME'])):
			#print("i=",i)
			if len(pdbstruct_new[chain]['ATOMNUM']) == 0:
				pdbstruct_new[chain]['ATOMNUM'].append(atomnum)
				pdbstruct_new[chain]['RESID'].append(resnum)
				atomnum += 1
				continue
			else:
				#print(pdbstruct[chain]['RESID'])
				#print(pdbstruct[chain]['RESID'][i-1])
				if pdbstruct[chain]['RESID'][i] == pdbstruct[chain]['RESID'][i-1]: # Check if the current res num is same as the prev res num
					pdbstruct_new[chain]['ATOMNUM'].append(atomnum)
					pdbstruct_new[chain]['RESID'].append(resnum)
					atomnum += 1
				else:
					resnum += 1
					pdbstruct_new[chain]['ATOMNUM'].append(atomnum)
					pdbstruct_new[chain]['RESID'].append(resnum)
					atomnum += 1
		resnum += 1			
	return pdbstruct_new				

