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
import re 
import itertools as it
import shutil
import glob
sys.stdout.flush()

"""
Parse the PDB file into a datastructure
"""
def parse_pdb(filename):
	# parse the contents of a PDB file as a dictionary. This is especially useful
	# if the PDB file has more than one chain.
	PDB_struct = {}
	fh = open(filename, 'r')
	for line in fh: 
		if re.match('^ATOM', line):
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
	
			res_name = line[17:20]
			res_id = line[22:26]
	
			res_name = re.sub('\s+','',res_name)
			res_id = re.sub('\s+','',res_id)
			res_id = int(res_id)
			atomname = line[13:15]
			if chain not in PDB_struct.keys():
				PDB_struct[chain] = dict()
				PDB_struct[chain]['COORD'] = []
				PDB_struct[chain]['RESNAME'] = []
				PDB_struct[chain]['RESID'] = []
				PDB_struct[chain]['ATOM_NAME'] = []
			PDB_struct[chain]['COORD'].append([X_coord, Y_coord, Z_coord])
			PDB_struct[chain]['RESNAME'].append(res_name)
			PDB_struct[chain]['RESID'].append(res_id)
			PDB_struct[chain]['ATOM_NAME'].append(atomname)
	fh.close()
	return PDB_struct

"""
Get the residue sequence for chain A and B (in this same order A and B) for
the given structure.
"""
def get_pdb_seq(pdb_dict, chain_order):
	seq = ''
	# Note 1W85 has all atoms missing for residue A-359, except for 'N'. Might need to double check the
	# results for this structure.
	for chain in chain_order:
		ca_pos = list(np.where(np.array(pdb_dict[chain]['ATOM_NAME']) == 'CA')[0])
		res_all = pdb_dict[chain]['RESNAME']
		res_ca = [res_all[i] for i in ca_pos] # Get the residue names corresponding only to the C-alpha atoms
		aa_seq = map_resname_to_id(res_ca) # Get the sequence of residues
		seq += aa_seq
	return seq

"""
Get the residue ids from PDB data structure
"""
def get_pdb_res_id(pdb_dict, chain_order):
	res_ids = []
	for chain in chain_order:
		ca_pos = list(np.where(np.array(pdb_dict[chain]['ATOM_NAME']) == 'CA')[0])
		res_ids_all = pdb_dict[chain]['RESID']
		res_ids_ca = [res_ids_all[i] for i in ca_pos] # Get the residue names corresponding only to the C-alpha atoms
		res_ids.extend(res_ids_ca)
	return res_ids	

"""
Convert an array of 3-lettered amino acid sequence to its single alphabet sequence
"""
def map_resname_to_id(res_name_seq):
	resname_2_id = { 
	    'ALA' : 'A',
	    'ARG' : 'R',
	    'ASN' : 'N',
	    'ASP' : 'D',
	    'CYS' : 'C',
	    'GLY' : 'G',
	    'GLN' : 'Q',
	    'GLU' : 'E',
	    'HIS' : 'H',
	    'ILE' : 'I',
	    'LEU' : 'L',
	    'LYS' : 'K',
	    'MET' : 'M',
	    'PRO' : 'P',
	    'PHE' : 'F',
	    'SER' : 'S',
	    'THR' : 'T',
	    'TRP' : 'W',
	    'TYR' : 'Y',
	    'VAL' : 'V',}
	res_seq = ''.join([resname_2_id[res_name] for res_name in res_name_seq]) # use list comprehension to map the 3-letter residue code to its alphabet
	return res_seq

"""
Get the native pdb sequence from the ccmpred alignment file
"""
def get_ccmpred_native_seq(ccmpred_aln_file):
	aln_seq = ''
	data = [line.rstrip() for line in open(ccmpred_aln_file)]
	aln_seq = data[0]
	return aln_seq	

"""
Align the ccmpred and pdb native sequences
"""
def align_seqs(ccmpred_seq, pdb_seq, aln_dir):
	# Remove any existing gaps from the ccmpred_seq
	ccmpred_seq = ccmpred_seq.replace('-','')
	ccmpred_seq = ccmpred_seq.replace(' ','')

	# Write these into a temporary alignment file and perform the alignment
	aln_input_file = aln_dir + '-' + 'ccmpred-pdb.aln'
	fh1 = open(aln_input_file, 'w')
	fh1.write(">aln_seq\n" + ccmpred_seq + "\n" + ">pdb_seq\n" + pdb_seq)
	fh1.close()
	aln_output_file = aln_dir + '-' +'ccmpred-pdb.aln.out'
	
	# Perform the alignment between the PDB sequence from the structure and the alignment file
	clust_o_path = os.environ['CLUSTALO']
	if not clust_o_path.endswith('/'):
		clust_o_path += '/'
	clust_o = clust_o_path + 'clustalo'	
	aln_cmd = clust_o + ' --force --outfmt=clu --wrap=100000 -i ' + aln_input_file + ' -o ' + aln_output_file # force overwrite any existing file and do not wrap
	os.system(aln_cmd)

	# Parse the alignment output file
	ccmpred_pdb_aln = [line.rstrip() for line in open(aln_output_file,'r')]
	return ccmpred_pdb_aln	

"""
Map the sequence in the PDB file to the sequence in the alignment file
A very strong underlying assumption in this function is that
all the residues in the pdb file are numbered serially. The residue numbers
must be unique, otherwise, this function will not work.
"""
def map_pdb_ccmpred_pos(ccmpred_pdb_aln, ccmpred_aln_seq_global, pdb_res_ids):
	ccmpred_aln = ''
	pdb_aln = ''
	for line in ccmpred_pdb_aln:
		if line.startswith('aln_seq'):
			split_line = line.split()
			ccmpred_aln = split_line[1]
		elif line.startswith('pdb_seq'):	
			split_line = line.split()
			pdb_aln = split_line[1]
	ccmpred_aln = ccmpred_aln.replace(' ','')
	pdb_aln = pdb_aln.replace(' ','')
	ccmpred_seq = ccmpred_aln.replace('-','')		
	pdb_seq = pdb_aln.replace('-','')
	mapping_dict1 = dict()
	mapping_dict2 = dict()
	pdb_seq_pos = -1
	ccmpred_seq_pos = -1
	### Very Tricky Business - Mapping the residue indices from the PDB sequence space to the CCMPred alignment space! ###
	# Stage 1 (Preliminary mapping): Map the positions in the pdb seq to the aligned ccmpred seq (with the pdb seq).
	for i in range(0,len(pdb_aln)):
		res_pdb = pdb_aln[i]
		res_ccmpred = ccmpred_aln[i]
		if res_pdb != '-' and res_ccmpred != '-': # Both positions have aa
			pdb_seq_pos += 1
			ccmpred_seq_pos += 1
			pdb_res = pdb_seq[pdb_seq_pos]
			ccmpred_res = ccmpred_seq[ccmpred_seq_pos]
			current_id = pdb_res_ids[pdb_seq_pos]
			mapping_dict1[current_id] = ccmpred_seq_pos
		elif res_pdb != '-' and res_ccmpred == '-': # residue in pdb was filtered in ccmpred
			pdb_seq_pos +=1 
			pdb_res = pdb_seq[pdb_seq_pos]
			current_id = pdb_res_ids[pdb_seq_pos]
			mapping_dict1[current_id] = ''
		elif res_pdb == '-' and res_ccmpred != '-': # residue in ccmpred and is missing in pdb
			ccmpred_seq_pos += 1
			ccmpred_res = ccmpred_seq[ccmpred_seq_pos]
	
	# Stage 2 (Final mapping): Map the positions in the pdb seq to the aligned ccmpred seq (with all the other sequences)
	# using results from stage 1.
	res_id = list(mapping_dict1.keys())
	res_id.sort() # sort in increasing order
	for r in res_id:
		if mapping_dict1[r] == '':
			continue
		else:
			# Find the corresponding position for aln_pos on ccmpred_aln_seq_global
			aln_pos = int(mapping_dict1[r]) # The residue index on the ccmpred sequence aligned to the PDB seq 
											# after excluding gaps 
			current_res_index = -1
			# We are trying to find the index of the residue on the ccmpred global alignment (alignment of the
			# ccmpred native seq with all the other sequences) for aln_pos
			# which corresponds to the 
			for n in range(0,len(ccmpred_aln_seq_global)):
				if ccmpred_aln_seq_global[n] != '-':
					current_res_index += 1
					if current_res_index == aln_pos: # This is the position corresponding to aln_pos and so, we record the n for this position
						mapping_dict2[r] = n
						break
	return mapping_dict2					
"""
Parse the kfc file and return the hotspot information
"""
def parse_kfc(hotspotfile,pred_class):
	hotspot_dict = {} # to store the parsed kfc information.
	flag = 0 
	fh = open(hotspotfile,'r')
	for line in fh: 
		if re.match('^.*?ChainSet1', line): # Capture the chain assignments
			chains = re.findall('^.*?ChainSet1:(.*?)\s+?ChainSet2:(.*?$)',line)
			set1 = chains[0][0]
			set2 = chains[0][1]
			set1 = set1.replace(" ","")
			set2 = set2.replace(" ","")
			hotspot_dict['CHAINSET1'] = set1
			hotspot_dict['CHAINSET2'] = set2
		if re.match('^\s*?-+?\s*$', line): # the hot spot information begins after a '----' line. Use this as a marker.
			flag = 1;
		elif flag == 1 and not re.match('^\s+?$', line): 
			line_elems = line.split();
			chain = line_elems[0]
			res = int(line_elems[2])
			if pred_class == 1: # Hotspot residues predicted by KFC2A
				category = line_elems[3] 
			elif pred_class == 2: # Hotspot residues predicted by KFC2B
				category = line_elems[5] 
			elif pred_class == 3: # Residue is hotspot if predicted as hotspot by both KFC2A and KFC2B
				if line_elems[3] == 'Hotspot' and line_elems[5] == 'Hotspot':
					category = 'Hotspot'
				else:
					category = ''
			else:		
				sys.exit("Error!!! Non-standard required argument type \'pred_class\' encountered!!!\n")
			if category == 'Hotspot':
				if chain not in hotspot_dict.keys(): 
					hotspot_dict[chain] = [] # intiate an empty list for a chain that is was previously not included in the dict
				hotspot_dict[chain].append(res)
		else:
			continue
	fh.close()  
	return hotspot_dict

"""
Parse the set of unique contacts
"""
def parse_contacts(contact_file):
	contact_dict = dict()
	fh = open(contact_file, 'r')
	for line in fh: 
	    if re.match(line, '^\s+?$'):
	        continue
	    else:
	        line_elems = line.split(',')
	        res1 = int(line_elems[0]) # residue 1 is always from chain A
	        res2 = int(line_elems[1]) # residue 2 is always from chain B
	        chain1 = 'A' # Based on our naming conventions, chain 1 is always chain A
	        chain2 = 'B' # Based on our naming conventions, chain 2 is always chain B
	        if chain1 not in contact_dict.keys():
	            contact_dict[chain1] = dict()
	        if chain2 not in contact_dict.keys():   
	            contact_dict[chain2] = dict()
	        if res1 not in contact_dict[chain1].keys():
	            contact_dict[chain1][res1] = []
	        if res2 not in contact_dict[chain2].keys(): 
	            contact_dict[chain2][res2] = []
	        contact_dict[chain1][res1].append(res2)
	        contact_dict[chain2][res2].append(res1)
	fh.close()
	return contact_dict

"""
Score a given binding mode using the KFC-E scoring function
"""
def score_hotspots(chain_id, hotspot_res, contact_dict, contact_map, mapped_pos_dict, cons_scores):
	total_score = 0
	rv = 0
	if chain_id in contact_dict.keys():
		for res_i in hotspot_res:
			contact_res_in_partner = [] # Store the contacting residues in the partner chain
			# Get the contacting residues for res_i
			if res_i in contact_dict[chain_id].keys():
				contact_res = contact_dict[chain_id][res_i]
				rv = 1 # In case we need to check if the hotspot residues were found in the contact dictionary
				partner_res_coevol_score = 0;
				if res_i in mapped_pos_dict.keys():
					pos_i = mapped_pos_dict[res_i] # Get the corresponding position for res_i in the ccmpred_alignment
					res_i_coevol_score = 0 # We will score each hotspot residue based on our conditional scoring function and 
									# then add it to the total score
					for res_j in contact_res:
						if res_j in mapped_pos_dict.keys(): # only if res_j is included in the alignment
							pos_j = mapped_pos_dict[res_j] # Get the corresponding position for res_j in the ccmpred_alignment
							score_ij = contact_map[pos_i, pos_j]
							res_i_coevol_score += score_ij
							contact_res_in_partner.append(res_j)
						else: 	# spit out a warning message and continue
							continue
					# If there are contacting residues in the interacting chain, then get the strength of co-evolution between them
					res_i_cons_score = cons_scores[pos_i] # Get the conservation of residue i
					if len(contact_res_in_partner) > 2:
						partner_res_comb = it.combinations(contact_res_in_partner,2)
						num_pairs = 0;
						for comb_i in partner_res_comb:
							r1 = mapped_pos_dict[comb_i[0]]
							r2 = mapped_pos_dict[comb_i[1]]
							num_pairs += 1
							partner_res_coevol_score += contact_map[r1,r2]
						if partner_res_coevol_score < 0.00000001: # Only for really wierd cases for which the contacts have very very low coevolution
							partner_res_coevol_score = 0.00000001
						conditional_score = (res_i_coevol_score * res_i_cons_score)/partner_res_coevol_score # conditional scheme 1 (results_3 dir) # works best
						total_score += conditional_score
					else:
						conditional_score = res_i_coevol_score * res_i_cons_score
						total_score += conditional_score
				else:
					continue
					#print ("Warning!!! Residue ", res_i, " not found in mapped pos dict!!! Most likely that the residue is there in the PDB but, omitted from the alignment!")
			else:
				continue
				#print ("No contact found for res ", res_i, " in chain ", chain_id)
	return (total_score, rv)

"""
Caller function to score binding modes with the KFC-E scoring function
"""
def scoreposes(nativepdbfile, chainorder, ccmpredalnfile, ccmpredoutfile, consfile, posedir, contactdir, hotspotdir, tmp_aln_dir, kfc_pred_class, outputdir):
	# Get the sequence for the native complex from the PDB and the CCM Pred alignment files
	native_dict = parse_pdb(nativepdbfile)
	pdb_seq = get_pdb_seq(native_dict, chainorder)
	pdb_res_ids = get_pdb_res_id(native_dict, chainorder)

	ccmpred_seq = get_ccmpred_native_seq(ccmpredalnfile) # The native aligned sequence (the first sequence) from the ccmpred alignment file
	nativepdb_prefix = nativepdbfile.split('.pdb')[0]
	# Align these sequences
	ccmpred_pdb_aln = align_seqs(ccmpred_seq, pdb_seq, tmp_aln_dir)
	#print ("ccmpred pdb aln : ", ccmpred_pdb_aln)
		
	# Get a mapping between the positions of the pdb seq and the alignment seq
	mapped_pos_dict = map_pdb_ccmpred_pos(ccmpred_pdb_aln, ccmpred_seq, pdb_res_ids)
		
	# For each sampled pose in the fraction contact file:
	#	1. Get the KFC2 hotspots
	#	2. Get the contacts for the hotspots
	#	3. Calculate cumulative hotspot coevolution scores
	chain_A_scores = []
	chain_B_scores = []
	chain_A_pose_ids = []
	chain_B_pose_ids = []
	contact_map = np.loadtxt(ccmpredoutfile, dtype=float)
	
	workdir = os.getcwd()
	os.chdir(posedir)
	posefiles = glob.glob('*.pdb')
	os.chdir(workdir)
	
	# Read the conservation scores into an array
	cons_scores = np.loadtxt(consfile, dtype=float)
	cons_scores = -cons_scores
	cons_scores = cons_scores + abs(min(cons_scores)) # rescale the conservation score to start from 0 and to ensure that
														# most conserved residues score high
	cons_scores = cons_scores / max(cons_scores) # re-scale the scores from 0 to 1													
	complex_name = ''
	for pose_i in posefiles:
		pose_id = int(pose_i.rsplit('_', 1)[1].split('.pdb')[0])
		pose_i_prefix = pose_i.split('.pdb')[0]
		if not complex_name:
			complex_name = pose_i.rsplit('_', 1)[0]
		
		# Parse the KFC2 output for this pose and get the hotspots
		kfc_file = hotspotdir  + str(pose_i_prefix) + '.kfc.results'
		hotspot_dict = parse_kfc(kfc_file, kfc_pred_class)
			
		# Parse the contact file and get the list of interface contacts
		# Instead of using the .fuc file we will use the .fc file as the .fuc file
		# is often missing some contacts for some complexes.
		contact_file = contactdir + pose_i_prefix + '.fuc'
		contact_dict = parse_contacts(contact_file)
			
		# Calculate cumulative coevolution scores for the hotspot residues
		chain1 = 'A' # Assume that the protein has only two chains - chain A and B
		chain2 = 'B'
		disp_str = "*** " + complex_name + " Pose : " + str(pose_id) +  "***\n"
	
		# We will only process those poses which have predicted hotspots on atleast one of the chains. For all other poses,
		# we will assign a score of 0.
		if chain1 not in hotspot_dict.keys() and chain2 not in hotspot_dict.keys():
			chain_A_scores.append(0)
			chain_B_scores.append(0)
			chain_A_pose_ids.append(pose_id)
			chain_B_pose_ids.append(pose_id)
			continue
		if chain1 in hotspot_dict.keys():
			chain1_hotspots = hotspot_dict[chain1]
			[chain_1_score, rv] = score_hotspots(chain1, chain1_hotspots, contact_dict, contact_map, mapped_pos_dict, cons_scores)
			chain_A_scores.append(chain_1_score)
			chain_A_pose_ids.append(pose_id)
			disp_str += ", Chain A score = " + str(round(chain_1_score,4))
		else:
			chain_A_scores.append(0)
			chain_A_pose_ids.append(pose_id)
			disp_str += ", Chain A score = " + str(0)
		if chain2 in hotspot_dict.keys():
			chain2_hotspots = hotspot_dict[chain2]
			[chain_2_score, rv] = score_hotspots(chain2, chain2_hotspots, contact_dict, contact_map, mapped_pos_dict, cons_scores)
			chain_B_scores.append(chain_2_score)
			chain_B_pose_ids.append(pose_id)
			disp_str += ", Chain B score = " + str(round(chain_2_score,4))
		else:
			chain_B_scores.append(0)
			chain_B_pose_ids.append(pose_id)
			disp_str += ", Chain B score = " + str(0)
	
	# Get the cumulative scores by adding the individual scores for chain A and B and then, normalizing
	chain_AB_scores = np.array(chain_A_scores) + np.array(chain_B_scores)
	chain_AB_scores = list(chain_AB_scores/max(chain_AB_scores)) # Store the cumulative scores for chains A and B

	# Write the scores into a score file
	score_file = outputdir + complex_name + '_KFC_E_scores.csv'
	matrix_AB = np.array([chain_A_pose_ids, chain_AB_scores]).transpose()

	# Sort by co-evolution scores
	matrix_AB_sorted = matrix_AB[matrix_AB[:,1].argsort()[::-1]]
	np.savetxt(score_file, matrix_AB_sorted, fmt='%0.4f')
	return score_file, matrix_AB

"""
Wrapper function for scoring the poses
"""
def scoreposes_wrapper(nativepdbfile, ccmpredalnfile, ccmpredoutfile, consfile, posedir, contactdir, hotspotdir, outputdir):
	chainorder = ['A', 'B']
	tmp_aln_dir = outputdir + 'tmp_aln'
	if os.path.isdir(tmp_aln_dir):
		shutil.rmtree(tmp_aln_dir)
	try:
		os.mkdir(tmp_aln_dir)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass
	kfc_pred_class = 1
	return scoreposes(nativepdbfile, chainorder, ccmpredalnfile, ccmpredoutfile, consfile, posedir, contactdir, hotspotdir, tmp_aln_dir, kfc_pred_class, outputdir)

