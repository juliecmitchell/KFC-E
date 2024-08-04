# KFC-E
Instructions for running KFC-E

The code for KFC-E has been tested on a machine with 64 bit Ubuntu (version 18.04) and is expected to be used
in a Linux machine with similar OS configuration.


A. Software
Download and install the following software.
	1. KFC2 (https://mitchell-lab.biochem.wisc.edu/KFC_Server/index.php)
	2. Rate4Site (https://m.tau.ac.il/~itaymay/cp/rate4site.html)
	3. Clustal Omega (http://www.clustal.org/omega/)
	4. ZDOCK (http://zdock.umassmed.edu/software/)
	5. CCMPRED (https://github.com/soedinglab/CCMpred)
	6. Python 3.6 (Preferably with the Anaconda package)
	7. Python 2.7

B. Environment variables
Set the following environment variables
	1. KFC2 (Set to the path having the KFC2 executable. The executable must be named as kfc2.sh.)
	2. CCMPRED (Set to the path having the ccmpred executable. The executable must be named 'ccmpred'.)
	3. RATE4SITE (Set to the path having the rate4site executable. The executable must be named 'rate4site'.)
	4. ZDOCK (Set to the path having the ZDOCK executable. The executable must be named 'zdock'.)
	5. CLUSTALO (Set to the path having the Clustal Omega executable. The executable must be named clustalo.)
	6. KFCE_HOME (Set to the KFCE folder)

C. Compile the contacts 'C' script and rename the executable to contacts

D. Required files
	a. Paired alignment file: This should be in the CCMpred format, other wise use the 'convert_alignment.py' script inside
	the CCMpred folder to convert.
	b. PDB file with 2 chains

E. Running KFC-E
	Example: python3.6 KFC_E.py --inputdir test_input --outputdir test_output --pdbfile 1rm6_BC.pdb --alnfile 1RM6_B_C_AB.fas
	Note: The pdb and the alignment files must be placed inside the test_input directory.

