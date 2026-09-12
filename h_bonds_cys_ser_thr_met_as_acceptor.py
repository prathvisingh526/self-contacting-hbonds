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

#^########################################################################################################################################################################################################################################################################
def calculator_and_writer(structure, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, relevant_neighbor_dist_cutoff, pdb_id_molprobity, is_cys_disulfide_bonded, molprobity_files, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, 
			  dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_donor_atom, desired_hydrogen_atom):
	distance_primary_sequence = desired_donor_atom.get_parent().get_full_id()[3][1]        - desired_acceptor_atom.get_parent().get_full_id()[3][1]
	distance_DA               = np.linalg.norm(np.array(desired_donor_atom.get_coord())    - np.array(desired_acceptor_atom.get_coord()))
	distance_HA               = np.linalg.norm(np.array(desired_hydrogen_atom.get_coord()) - np.array(desired_acceptor_atom.get_coord()))
	angle_DHA                 = math.degrees(calc_angle(desired_donor_atom.get_vector(),    desired_hydrogen_atom.get_vector(), desired_acceptor_atom.get_vector()))
	angle_HAAn                = math.degrees(calc_angle(desired_hydrogen_atom.get_vector(), desired_acceptor_atom.get_vector(), desired_antecedent_atom1.get_vector()))

	if desired_donor_atom.get_name() in dict_donor_and_hydrogen_atoms_ON.keys() and desired_acceptor_atom.get_name() == "SD" and desired_acceptor_atom.get_parent().get_resname() == "MET":
		vector_accep_hydrogen = np.array(desired_hydrogen_atom.get_coord())    - np.array(desired_acceptor_atom.get_coord()) #vector SD-H
		vector_accep_ant1     = np.array(desired_antecedent_atom1.get_coord()) - np.array(desired_acceptor_atom.get_coord()) #vector SD-CG
		vector_accep_ant2     = np.array(desired_antecedent_atom2.get_coord()) - np.array(desired_acceptor_atom.get_coord()) #vector SD-CE
		vector_bisector       = vector_accep_ant1 + vector_accep_ant2                                                        #vector which is bisector of vector SD-CG and vector SD-CE

		vector_accep_hydrogen_normalized = vector_accep_hydrogen / np.linalg.norm(vector_accep_hydrogen)
		vector_bisector_normalized       = vector_bisector       / np.linalg.norm(vector_bisector)

		dot_product                                = np.dot(vector_accep_hydrogen_normalized, vector_bisector_normalized) #dot product of vector SD-H and bisector vector (required to calculate the angle between them)
		magnitude_vector_accep_hydrogen_normalized = np.linalg.norm(vector_accep_hydrogen_normalized)                     #magnitude of vector SD-H
		magnitude_vector_bisector_normalized       = np.linalg.norm(vector_bisector_normalized)                           #magnitude of bisector vector
		angle_rad                                  = np.arccos(dot_product / (magnitude_vector_accep_hydrogen_normalized * magnitude_vector_bisector_normalized)) #angle between vector SD-H and bisector vector
		angle_deg                                  = np.degrees(angle_rad)
		angle_HAAn                                 = angle_deg

	if segments_relevant_chain != None: #! it is equal to "None" when running for globular dataset
		flag2 = False #flag1 has been used in the acceptor code which sends data here
		for segment_range in segments_relevant_chain["HBCIL"]:
			lower_limit_of_segment_range = int(segment_range.split("-")[0])
			upper_limit_of_segment_range = int(segment_range.split("-")[1])
			if desired_donor_atom.get_parent().get_full_id()[3][1] in range(lower_limit_of_segment_range, upper_limit_of_segment_range + 1):
				flag2 = True

	if desired_donor_atom.get_name() in dict_donor_and_hydrogen_atoms_ON.keys() and desired_acceptor_atom.get_name() in dict_acceptor_and_antecedent_atoms_ON.keys():
		distance_DA_geom = dict_geom_crit[donor_atoms_ON][acceptor_atoms_ON][0]
		distance_HA_geom = dict_geom_crit[donor_atoms_ON][acceptor_atoms_ON][1]
		angle_DHA_geom   = dict_geom_crit[donor_atoms_ON][acceptor_atoms_ON][2]
		angle_HAAn_geom  = dict_geom_crit[donor_atoms_ON][acceptor_atoms_ON][3]

	if desired_donor_atom.get_name() in dict_donor_and_hydrogen_atoms_ON.keys() and desired_acceptor_atom.get_name() in dict_acceptor_and_antecedent_atoms_S.keys():
		distance_DA_geom = dict_geom_crit[donor_atoms_ON][acceptor_atoms_S][0]
		distance_HA_geom = dict_geom_crit[donor_atoms_ON][acceptor_atoms_S][1]
		angle_DHA_geom   = dict_geom_crit[donor_atoms_ON][acceptor_atoms_S][2]
		angle_HAAn_geom  = dict_geom_crit[donor_atoms_ON][acceptor_atoms_S][3]

	if desired_donor_atom.get_name() in dict_donor_and_hydrogen_atoms_S.keys() and desired_acceptor_atom.get_name() in dict_acceptor_and_antecedent_atoms_ON.keys():
		distance_DA_geom = dict_geom_crit[donor_atoms_S][acceptor_atoms_ON][0]
		distance_HA_geom = dict_geom_crit[donor_atoms_S][acceptor_atoms_ON][1]
		angle_DHA_geom   = dict_geom_crit[donor_atoms_S][acceptor_atoms_ON][2]
		angle_HAAn_geom  = dict_geom_crit[donor_atoms_S][acceptor_atoms_ON][3]

	if desired_donor_atom.get_name() in dict_donor_and_hydrogen_atoms_S.keys() and desired_acceptor_atom.get_name() in dict_acceptor_and_antecedent_atoms_S.keys():
		distance_DA_geom = dict_geom_crit[donor_atoms_S][acceptor_atoms_S][0]
		distance_HA_geom = dict_geom_crit[donor_atoms_S][acceptor_atoms_S][1]
		angle_DHA_geom   = dict_geom_crit[donor_atoms_S][acceptor_atoms_S][2]
		angle_HAAn_geom  = dict_geom_crit[donor_atoms_S][acceptor_atoms_S][3]

	if distance_DA <= distance_DA_geom and distance_HA <= distance_HA_geom and angle_DHA >= angle_DHA_geom and angle_HAAn >= angle_HAAn_geom:
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor atom",                                         []).append(desired_donor_atom.get_name())                          #donor's atom name
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor resname",                                      []).append(desired_donor_atom.get_parent().get_resname())          #donor's residue name
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor resno",                                        []).append(desired_donor_atom.get_parent().get_full_id()[3][1])    #donor's residue no
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor reschain",                                     []).append(desired_donor_atom.get_parent().get_full_id()[2])       #donor's parent chain
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor PDB ID",                                       []).append(pdb_id_molprobity[0:4].upper())                         #donor's parent molecule
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor insertion code",                               []).append(desired_donor_atom.get_parent().get_full_id()[3][2])    #donor's insertion code

		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor atom",                                      []).append(desired_acceptor_atom.get_name())                       #acceptor's atom name
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor resname",                                   []).append(desired_acceptor_atom.get_parent().get_resname())       #acceptor's residue name
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor resno",                                     []).append(desired_acceptor_atom.get_parent().get_full_id()[3][1]) #acceptor's residue no
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor reschain",                                  []).append(desired_acceptor_atom.get_parent().get_full_id()[2])    #acceptor's parent chain
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor PDB ID",                                    []).append(pdb_id_molprobity[0:4].upper())                         #acceptor's parent molecule
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor insertion code",                            []).append(desired_acceptor_atom.get_parent().get_full_id()[3][2]) #acceptor's insertion code
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor CYS disulfide bonded",                      []).append(is_cys_disulfide_bonded)

		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("H-atom involved",                                    []).append(desired_hydrogen_atom.get_name())                       #if the donor atom is bound to two H-atoms, this will help you know which among the two is engaged in H-bonding
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("DA distance (Å)",                                    []).append(round(float(distance_DA), 2))
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("HA distance (Å)",                                    []).append(round(float(distance_HA), 2))
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("DHA angle (degrees)",                                []).append(round(angle_DHA, 2))
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("HAAn angle/HABisector angle (degrees)",              []).append(round(angle_HAAn, 2))
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Occupancy of donor atom",                            []).append(desired_donor_atom.get_occupancy())
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Occupancy of acceptor atom",                         []).append(desired_acceptor_atom.get_occupancy())
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Occupancy of hydrogen atom",                         []).append(desired_hydrogen_atom.get_occupancy())
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("DA distance in primary sequence",                    []).append(distance_primary_sequence)
		dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Resolution (Å)",                                     []).append(dict_resolutions[pdb_id_molprobity[0:4].upper()])
		if segments_relevant_chain != None:
			dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor residue belongs to TM segment of the same chain", []).append(flag2)

		#^#######################################################################################################################
		#^ SSEC code (identical for "tm" and "globular" datasets):
		#^#######################################################################################################################
		stride_return_value = ssec_strc_phi_psi_chi1_sasa_value_returner.sec_str_value_returner(desired_donor_atom, desired_acceptor_atom, pdb_id_molprobity, molprobity_files)
		if stride_return_value == False: # this happens when the stride file corresponding to the PDB structure isn't found
			dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor sec. str.",    []).append("Stride file not found") # I have used "Stride file not found" and not "NA" for a reason
			dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor sec. str.", []).append("Stride file not found")
		if stride_return_value != False:
			dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Donor sec. str.",    []).append(stride_return_value["Donor ssec"])
			dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault("Acceptor sec. str.", []).append(stride_return_value["Acceptor ssec"])

		#^###############################################################################################################################
		#^ PHI PSI CHI1 code (identical for "tm" and "globular" datasets):
		#^###############################################################################################################################
		phi_psi_chi1_return_value = ssec_strc_phi_psi_chi1_sasa_value_returner.chimera_phi_psi_chi1_value_returner(desired_donor_atom, desired_acceptor_atom, pdb_id_molprobity, molprobity_files)
		for key, value in phi_psi_chi1_return_value.items(): # it will look like: {"phi" : {"donor" : "None", "acceptor" : "None"}, {"psi" : {"donor" : "None", "acceptor" : "None"}, {"chi1" : {"donor" : "None", "acceptor" : "None"}}
			for acc_don, anglez in value.items():
				dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault(key + " angle " + acc_don + " (degrees)",    []).append(anglez)

		#^##################################################################
		#^ SASA code (the following code is ONLY for "globular" dataset):
		#^##################################################################
		sasa_return_value = ssec_strc_phi_psi_chi1_sasa_value_returner.sasa_value_returner(desired_donor_atom, desired_acceptor_atom, pdb_id_molprobity, molprobity_files)
		for key, value in sasa_return_value.items():
			dictdata[desired_acceptor_atom.get_parent().get_resname()].setdefault(key, []).append(value)

#^########################################################################################################################################################################################################################################################################
def h_bonds_finder(structure, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, relevant_neighbor_dist_cutoff, pdb_id_molprobity, is_cys_disulfide_bonded, molprobity_files, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, 
		   dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain):
	atoms     = list(structure.get_atoms())    #I had to convert the generator to list because iterating over generator throws error
	neighbors = list()
	neighbors = NeighborSearch(atoms).search(desired_acceptor_atom.coord, relevant_neighbor_dist_cutoff)
	if len(neighbors) == 0:
		dict_no_neighbors.setdefault(pdb_id_molprobity[0:4], desired_acceptor_atom.get_parent().get_resname() + str(desired_acceptor_atom.get_parent().get_full_id()[3][1]) + desired_acceptor_atom.get_parent().get_full_id()[2])
	if len(neighbors) != 0:
		for neighbor in neighbors:
			if (neighbor.get_parent().get_full_id()[3][0] == " " and # if neighbor atom's parent's residue is a non-heteroatom entity 
				neighbor.get_parent().get_full_id()[3][2] == " " and # if neighbor atom's parent's residue has no insertion code
				neighbor.get_parent().get_full_id()[2] == desired_acceptor_atom.get_parent().get_full_id()[2]): #if neighbor atom's parent's residue belongs to the same chain as the donor atom

				merged_dict = {**dict_donor_and_hydrogen_atoms_ON, **dict_donor_and_hydrogen_atoms_S}
				if neighbor.get_name() in merged_dict.keys():
					desired_donor_atom = neighbor
					check_flag         = False

					if desired_acceptor_atom.get_name() not in merged_dict.keys(): # if the acceptor atom is not in "merged_dict.keys()", it means that it isnt bonded to a H-atom (this condition will will be true for MET but not for CYS SER THR)
						check_flag = True

					if desired_acceptor_atom.get_name() in merged_dict.keys():     # if the acceptor is in "merged_dict.keys()", grab all the H-atoms bonded to it
						h_atoms_bonded_to_acceptor_atom = list()
						for check_atom in desired_acceptor_atom.get_parent().get_atoms():
							if check_atom.get_name() in merged_dict[desired_acceptor_atom.get_name()]: # grab all the H-atoms bonded to the acceptor atom
								h_atoms_bonded_to_acceptor_atom.append(check_atom)

						if len(h_atoms_bonded_to_acceptor_atom) == 0: # if the acceptor atom should have H-atoms bonded to it but in reality there aren't any (e.g., CYS side chain SG atom)
							check_flag = True

						if len(h_atoms_bonded_to_acceptor_atom) != 0:
							# print(desired_acceptor_atom.get_name(), desired_acceptor_atom.get_parent().get_resname(), h_atoms_bonded_to_acceptor_atom)
							dict_check_distances = dict()
							dict_check_distances.setdefault("DA distance",  np.linalg.norm(np.array(desired_donor_atom.get_coord())  - np.array(desired_acceptor_atom.get_coord())))
							dict_check_distances.setdefault("DH distances", [np.linalg.norm(np.array(desired_donor_atom.get_coord()) - np.array(hatom.get_coord())) for hatom in h_atoms_bonded_to_acceptor_atom])
							if dict_check_distances["DA distance"] < min(dict_check_distances["DH distances"]):
								check_flag = True

					if check_flag == True: # only proceed either if the acceptor atom isn't bonded to a H-atom or if it is bonded to a H-atom but the DA distance is the least
						for don_h_atom in desired_donor_atom.get_parent().get_atoms():
							if don_h_atom.get_name() in merged_dict[desired_donor_atom.get_name()]:
								desired_hydrogen_atom = don_h_atom
								calculator_and_writer(structure, desired_acceptor_atom, desired_antecedent_atom1, desired_antecedent_atom2, relevant_neighbor_dist_cutoff, pdb_id_molprobity, is_cys_disulfide_bonded, molprobity_files, dict_geom_crit, dict_no_neighbors, dictdata, dict_resolutions, 
			      dict_donor_and_hydrogen_atoms_ON, dict_donor_and_hydrogen_atoms_S, dict_acceptor_and_antecedent_atoms_ON, dict_acceptor_and_antecedent_atoms_S, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S, segments_relevant_chain, desired_donor_atom, desired_hydrogen_atom)
