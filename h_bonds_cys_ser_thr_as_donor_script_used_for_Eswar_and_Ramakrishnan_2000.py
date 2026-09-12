import re
import os
import pprint
import json
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.PDB import *
import importlib
import sys
sys.path.append("G:\\My Drive\\python scripts and parameters\\")
import ssec_strc_phi_psi_chi1_sasa_value_returner
importlib.reload(ssec_strc_phi_psi_chi1_sasa_value_returner)

#########################################################################################################################################################################################################################################################################
def calculator_and_writer(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, dict_acceptor_and_antecedent_atoms_ON, 
				   dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, is_cys_disulfide_bonded):
	distance_primary_sequence = desired_donor_atom.get_parent().get_full_id()[3][1]        - desired_acceptor_atom.get_parent().get_full_id()[3][1]
	distance_DA               = np.linalg.norm(np.array(desired_donor_atom.get_coord())    - np.array(desired_acceptor_atom.get_coord()))
	distance_HA               = np.linalg.norm(np.array(desired_hydrogen_atom.get_coord()) - np.array(desired_acceptor_atom.get_coord()))
	angle_HDA                 = math.degrees(calc_angle(desired_hydrogen_atom.get_vector(), desired_donor_atom.get_vector(), desired_acceptor_atom.get_vector()))

	if distance_DA >= 2.4 and distance_DA <= 3.5 and angle_HDA >= 70 and angle_HDA <= 150:
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Donor atom",                            []).append(desired_donor_atom.get_name())                          #donor's atom name
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Donor resname",                         []).append(desired_donor_atom.get_parent().get_resname())          #donor's residue name
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Donor resno",                           []).append(desired_donor_atom.get_parent().get_full_id()[3][1])    #donor's residue no
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Donor reschain",                        []).append(desired_donor_atom.get_parent().get_full_id()[2])       #donor's parent chain
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Donor PDB ID",                          []).append(pdb_id_molprobity[0:4].upper())                         #donor's parent molecule
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Donor insertion code",                  []).append(desired_donor_atom.get_parent().get_full_id()[3][2])    #donor's insertion code

		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor atom",                         []).append(desired_acceptor_atom.get_name())                       #acceptor's atom name
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor resname",                      []).append(desired_acceptor_atom.get_parent().get_resname())       #acceptor's residue name
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor resno",                        []).append(desired_acceptor_atom.get_parent().get_full_id()[3][1]) #acceptor's residue no
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor reschain",                     []).append(desired_acceptor_atom.get_parent().get_full_id()[2])    #acceptor's parent chain
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor PDB ID",                       []).append(pdb_id_molprobity[0:4].upper())                         #acceptor's parent molecule
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor insertion code",               []).append(desired_acceptor_atom.get_parent().get_full_id()[3][2]) #acceptor's insertion code
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Acceptor CYS disulfide bonded",         []).append(is_cys_disulfide_bonded)

		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("H-atom involved",                       []).append(desired_hydrogen_atom.get_name())                       #if the donor atom is bound to two H-atoms, this will help you know which among the two is engaged in H-bonding
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("DA distance (Å)",                       []).append(round(float(distance_DA), 2))
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("HA distance (Å)",                       []).append(round(float(distance_HA), 2))
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("HDA angle (degrees)",                   []).append(round(float(angle_HDA), 2)) 
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Occupancy of donor atom",               []).append(desired_donor_atom.get_occupancy())
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Occupancy of acceptor atom",            []).append(desired_acceptor_atom.get_occupancy())
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Occupancy of hydrogen atom",            []).append(desired_hydrogen_atom.get_occupancy())
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("DA distance in primary sequence",       []).append(distance_primary_sequence)
		dictdata[desired_donor_atom.get_parent().get_resname()].setdefault("Resolution (Å)",                        []).append(dict_resolutions[pdb_id_molprobity[0:4].upper()])

#########################################################################################################################################################################################################################################################################
def h_bonds_finder(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, dict_acceptor_and_antecedent_atoms_ON, 
		   dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain):
	atoms     = list(structure.get_atoms()) #I had to convert the generator to list because iterating over generator throws error
	neighbors = list()
	neighbors = NeighborSearch(atoms).search(desired_donor_atom.coord, relevant_neighbor_dist_cutoff)
	if len(neighbors) == 0:
		dict_no_neighbors.setdefault(pdb_id_molprobity[0:4], desired_donor_atom.get_parent().get_resname() + str(desired_donor_atom.get_parent().get_full_id()[3][1]) + desired_donor_atom.get_parent().get_full_id()[2])
	if len(neighbors) != 0:
		for neighbor in neighbors:
			if (neighbor.get_parent().get_full_id()[3][0] == " " and # if neighbor atom's parent's residue is a non-heteroatom entity 
       			neighbor.get_parent().get_full_id()[3][2] == " " and # if neighbor atom's parent's residue has no insertion code
				neighbor.get_parent().get_full_id()[2]    == desired_donor_atom.get_parent().get_full_id()[2]): # if neighbor atom's parent's residue belongs to the same chain as the donor atom

				if neighbor.get_name() in dict_acceptor_and_antecedent_atoms_ON.keys() and neighbor != desired_donor_atom:
					desired_acceptor_atom   = neighbor
					is_cys_disulfide_bonded = "NA"
					check_flag              = False
					check_flag              = min_dist_checker(dict_donor_and_hydrogen_atoms_ON, desired_donor_atom, desired_acceptor_atom)

					if check_flag == True: # it will be True either if the acceptor atom isn't bonded to a H-atom or if it is bonded to a H-atom but the DA distance is the least

						if desired_acceptor_atom.get_name() == "NE2" and desired_acceptor_atom.get_parent().get_resname() == "GLN":
							desired_antecedent_atom1 = next((acc_atom for acc_atom in desired_acceptor_atom.get_parent().get_atoms() if acc_atom.get_name() == dict_acceptor_and_antecedent_atoms_ON[desired_acceptor_atom.get_name()][0]), None)
							desired_antecedent_atom2 = "empty"
							if desired_antecedent_atom1 is not None:
								calculator_and_writer(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, 
			      dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, is_cys_disulfide_bonded)

						if desired_acceptor_atom.get_name() == "NE2" and desired_acceptor_atom.get_parent().get_resname() == "HIS":
							desired_antecedent_atom1 = next((acc_atom for acc_atom in desired_acceptor_atom.get_parent().get_atoms() if acc_atom.get_name() == dict_acceptor_and_antecedent_atoms_ON[desired_acceptor_atom.get_name()][1]), None)
							desired_antecedent_atom2 = "empty"
							if desired_antecedent_atom1 is not None:
								calculator_and_writer(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, 
					dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, is_cys_disulfide_bonded)

						if desired_acceptor_atom.get_name() != "NE2":
							desired_antecedent_atom1 = next((acc_atom for acc_atom in desired_acceptor_atom.get_parent().get_atoms() if acc_atom.get_name() == dict_acceptor_and_antecedent_atoms_ON[desired_acceptor_atom.get_name()][0]), None)
							desired_antecedent_atom2 = "empty"
							if desired_antecedent_atom1 is not None:
								calculator_and_writer(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, 
					dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, is_cys_disulfide_bonded)

				if neighbor.get_name() in dict_acceptor_and_antecedent_atoms_S.keys():
					desired_acceptor_atom   = neighbor
					is_cys_disulfide_bonded = "NA"
					check_flag              = False
					check_flag              = min_dist_checker(dict_donor_and_hydrogen_atoms_S, desired_donor_atom, desired_acceptor_atom)

					if check_flag == True:

						if desired_acceptor_atom.get_name() == "SD" and desired_acceptor_atom.get_parent().get_resname() == "MET":
							desired_antecedent_atom1 = next((acc_atom for acc_atom in desired_acceptor_atom.get_parent().get_atoms() if acc_atom.get_name() == dict_acceptor_and_antecedent_atoms_S[desired_acceptor_atom.get_name()][0]), None)
							desired_antecedent_atom2 = next((acc_atom for acc_atom in desired_acceptor_atom.get_parent().get_atoms() if acc_atom.get_name() == dict_acceptor_and_antecedent_atoms_S[desired_acceptor_atom.get_name()][1]), None)
							if desired_antecedent_atom1 is not None and desired_antecedent_atom2 is not None:
								calculator_and_writer(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, 
					dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, is_cys_disulfide_bonded)

						if desired_acceptor_atom.get_name() == "SG" and desired_acceptor_atom.get_parent().get_resname() == "CYS":
							desired_antecedent_atom1 = next((acc_atom for acc_atom in desired_acceptor_atom.get_parent().get_atoms() if acc_atom.get_name() == dict_acceptor_and_antecedent_atoms_S[desired_acceptor_atom.get_name()][0]), None)
							desired_antecedent_atom2 = "empty"
							is_cys_disulfide_bonded  = False
							reschain_plus_resno      = str(desired_acceptor_atom.get_parent().get_full_id()[2]) + "." + str(desired_acceptor_atom.get_parent().get_full_id()[3][1]) # "reschain" plus "resno"
							if reschain_plus_resno in list_disulfides:
								is_cys_disulfide_bonded = True
							if desired_antecedent_atom1 is not None:
								calculator_and_writer(structure, desired_donor_atom, desired_hydrogen_atom, relevant_neighbor_dist_cutoff, pdb_id_molprobity, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, list_disulfides, molprobity_files, dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, 
					dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, is_cys_disulfide_bonded)

#########################################################################################################################################################################################################################################################################
def min_dist_checker(dictionary, desired_donor_atom, desired_acceptor_atom):
	if desired_acceptor_atom.get_name() not in dictionary.keys(): # if the acceptor atom is in "dict_acceptor_and_antecedent_atoms_ON.keys()" but not in "dict_donor_and_hydrogen_atoms_ON.keys()", it means that it isnt bonded to a H-atom
		return True
	
	if desired_acceptor_atom.get_name() in dictionary.keys():     # if the acceptor atom is bound to a H-atom, then grab all the H-atoms bonded to it
		list_of_h_atoms = list()
		for check_atom in desired_acceptor_atom.get_parent().get_atoms():
			if check_atom.get_name() in dictionary[desired_acceptor_atom.get_name()]: # grab all the H-atoms bonded to the acceptor atom
				list_of_h_atoms.append(check_atom)

		if len(list_of_h_atoms) == 0: # if the acceptor atom should have H-atoms bonded to it but in reality there aren't any (e.g., CYS side chain SG atom)
			return True

		if len(list_of_h_atoms) != 0:
			# print(desired_acceptor_atom.get_name(), desired_acceptor_atom.get_parent().get_resname(), list_of_h_atoms)
			dict_check_distances = dict()
			dict_check_distances.setdefault("DA distance",  np.linalg.norm(np.array(desired_donor_atom.get_coord())  - np.array(desired_acceptor_atom.get_coord())))
			dict_check_distances.setdefault("DH distances", [np.linalg.norm(np.array(desired_donor_atom.get_coord()) - np.array(hatom.get_coord())) for hatom in list_of_h_atoms])
			if dict_check_distances["DA distance"] < min(dict_check_distances["DH distances"]): #if this condition isn't satisfied, the function will return None
				return True
