#NOTE: This script will be execute as part of the code which was meant to create systems for interaction energy calculations for cases where the acceptor atom is mainchain N-atom of either self/nonself residue.

import os
import shutil
from chimera import runCommand

# don_atoms         = "SG,HG,HG1"      ####################### mentioning the exact the names of the donor half atoms is very important
don_atoms         = "SG,HG,CB,HB1,HB2,HB3" ####################### mentioning the exact the names of the donor half atoms is very important 
path_hadded_files = "D:\\hadded_files\\CYS\\"
path_system_files = "D:\\system_files\\"
if not os.path.isdir(path_system_files):
	os.makedirs(path_system_files)

for haddedfile in os.listdir(path_hadded_files):
	pdb_id      = haddedfile.split("_")[0]
	don_atom    = haddedfile.split("_")[1]
	don_resname = haddedfile.split("_")[2]
	don_resno   = int(haddedfile.split("_")[3])
	chain_id    = str(haddedfile.split("_")[4])

	acc_atom    = haddedfile.split("_")[6]
	acc_resname = haddedfile.split("_")[7]
	acc_resno   = int(haddedfile.split("_")[8])

	don_resno   = str(don_resno)
	acc_resno   = str(acc_resno)

	don_chain_id_and_resno = don_resno  + "." + chain_id

	don_filename    = pdb_id + "_" + don_atom + "_" + don_resname + "_" + don_resno + "_" + chain_id + "_don.pdb"
	acc_filename    = pdb_id + "_" + acc_atom + "_" + acc_resname + "_" + acc_resno + "_" + chain_id + "_accep.pdb"
	donacc_filename = pdb_id + "_" + don_atom + "_" + don_resname + "_" + don_resno + "_" + chain_id + "+_" + acc_atom + "_" + acc_resname + "_" + acc_resno + "_" + chain_id + "_donaccep.pdb"

	###############################################################################################
	# Create the "D" system by keeping only the "D" system atoms:
	###############################################################################################
	runCommand("open " + path_hadded_files + haddedfile)
	runCommand("select:" + don_chain_id_and_resno + "@" + don_atoms)
	runCommand("select invert sel")
	runCommand("delete sel")
	runCommand("write 0 " + path_system_files + don_filename)
	runCommand("close 0")

	###############################################################################################
	# Create the "A" system by deleting the "D" system atoms:
	###############################################################################################
	runCommand("open " + path_hadded_files + haddedfile)
	runCommand("select:" + don_chain_id_and_resno + "@" + don_atoms)
	runCommand("delete sel")
	runCommand("write 0 " + path_system_files + acc_filename)
	runCommand("close 0")

	###############################################################################################
	# Copy the "D+A" system to the same location as "D" and "A" systems:
	###############################################################################################
	shutil.copy(path_hadded_files + haddedfile, path_system_files + donacc_filename)

runCommand("stop")