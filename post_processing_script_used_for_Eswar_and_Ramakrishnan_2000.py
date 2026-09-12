import pandas as pd
import re

def undesired_bifurcates_remover(df1, excel_filename, reskey, donor_atoms_ON, donor_atoms_S, acceptor_atoms_ON, acceptor_atoms_S):
    df1.sort_values("Donor PDB ID", inplace = True)
    df2                       = df1.copy() #this line causes the code to throw the warning: "A value is trying to be set on a copy of a slice from a DataFrame". It doesn't mean that something is wrong with your code. It simply means that any changes made in "df2" will not affect "df1"
    dict_common_acceptor_info = dict()
    rows_indices_unwanted_bif = list()
    for i in df1.itertuples():
        acceptor_atom     = df1.at[i[0], "Acceptor atom"]
        acceptor_resname  = df1.at[i[0], "Acceptor resname"]
        acceptor_resno    = df1.at[i[0], "Acceptor resno"]
        acceptor_reschain = df1.at[i[0], "Acceptor reschain"]
        acceptor_pdb_id   = df1.at[i[0], "Acceptor PDB ID"]
        tup1              = tuple() + (acceptor_atom, acceptor_resname, acceptor_resno, acceptor_reschain, acceptor_pdb_id) #using tuple and not list because a list cannot become key of a dictionary
        dict_common_acceptor_info.setdefault(tup1, []).append(i[0])                                                         #it will look like: {('O', 'ALA', 52, 'A', '3DB2') [1299, 1300]} for an acceptor which is mentioned in >1 row

    for key, value in dict_common_acceptor_info.items():
        if len(value) > 1: #if len(value) > 1, it means that there are >1 rows with the same values of acceptor info.
            dict_donor = dict()
            for i in df1.itertuples():
                if i[0] in value:
                    donor_resno    = df1.at[i[0], "Donor resno"]
                    donor_reschain = df1.at[i[0], "Donor reschain"]
                    donor_pdb_id   = df1.at[i[0], "Donor PDB ID"]
                    ha_distance    = df1.at[i[0], "HA distance (Å)"]
                    tup2           = tuple() + (donor_resno, donor_reschain, donor_pdb_id)
                    dict_donor.setdefault(tup2, []).append(ha_distance) #if two donors atoms are donating to the same residue, they can either belong to the same residue or different residues. If they belong to same residue, then their HA distances will append to the same key

            for key2, value2 in dict_donor.items():
                if len(value2) > 1:                    #if "len(value2) > 1", it means that >1 donor atoms, belonging to the same residue, are donating to the same acceptor atom
                    for hadistance in value2:          #in such a case, only keep the instance which is closer to the acceptor atom
                        if hadistance != min(value2):  #if the value in the current iteration is not the smallest value (I am doing this because I want the index number of those rows whose HA distance isn't minimum so that I can delete them later)
                            for j in df2.itertuples(): #then search for it corresponding row in the dataframe & grab its row number
                                donor_resno    = df2.at[j[0], "Donor resno"]
                                donor_reschain = df2.at[j[0], "Donor reschain"]
                                donor_pdb_id   = df2.at[j[0], "Donor PDB ID"]
                                ha_distance    = df2.at[j[0], "HA distance (Å)"]
                                if donor_resno == key2[0] and donor_reschain == key2[1] and donor_pdb_id == key2[2] and ha_distance == hadistance:
                                    rows_indices_unwanted_bif.append(j[0])

    if len(rows_indices_unwanted_bif) == 0:
        print("No unwanted instances of bifurcated H-bonds found.")
        df3 = pd.DataFrame() #empty dataframe

    if len(rows_indices_unwanted_bif) > 0:
        print("Unwanted instances of bifurcated H-bonds found. Their info will be written in a separate tab of the excel file and then their info from the main dataframe will be deleted.")
        df3 = df2.loc[rows_indices_unwanted_bif]
        df3.sort_values("Donor PDB ID", inplace = True)

        for i in rows_indices_unwanted_bif[::-1]: #this loop will do the same thing as "for i in range(len(df1)-1, -1, -1)"; the loop being used right now requires one step less
            df1 = df1.drop([i], axis = 0)
        df1 = df1.reset_index(drop = True)

##################################################################################################################
#Create a datafame where DA occupancy is 1:
##################################################################################################################

    df4 = df1[(df1["Occupancy of donor atom"] == 1) & (df1["Occupancy of acceptor atom"] == 1)]
    df4.sort_values("Donor PDB ID", inplace = True)

##################################################################################################################
#Create a datafame for self contacts:
##################################################################################################################

    df5 = df4[(df4["Donor resname"] == df4["Acceptor resname"]) & (df4["Donor resno"] == df4["Acceptor resno"]) & (df4["Donor reschain"] == df4["Acceptor reschain"])]
    df5 = df5.reset_index(drop = True) #if you don't reset index, then if a self contact was in the 7th row of df4 and it gets added to the 1st row of df5, its index will still stay 7

##################################################################################################################
#Write dataframes to excel file:
##################################################################################################################

    if "sidechain is the donor" in excel_filename:
        filter_word2 = " don. "
    if "sidechain is the acceptor" in excel_filename:
        filter_word2 = " accep. "

    list_of_worksheets = [
                            [df1, "All inst. undesired bif excl."], 
                            [df2,  "All inst. undesired bif incl."], 
                            [df3,  "Unwanted bif rows"], 
                            [df4,  "Instances DA occup one"], 
                            [df5,  "Self contact instances"], 
                        ]

    with pd.ExcelWriter("D:\\" + excel_filename, engine = "xlsxwriter") as writer:
        for df_name, worksheet_name in list_of_worksheets:
            df_name.to_excel(writer, sheet_name = worksheet_name, index = False, na_rep = "NA")
