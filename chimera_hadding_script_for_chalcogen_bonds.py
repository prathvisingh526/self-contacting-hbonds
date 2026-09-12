#NOTE: After this script has been run, you will observe that the "TER" row gets populated by the entry of the residue in the last line. This is completely normal & doesn't interferes with any subsequent step.
#NOTE: In cases where the acceptor atom is mainchain oxygen atom, I had to rename the acceptor residue to a non-conventional amino acid name & its C atom name to CA because otherwise chimera recognizes it as a 
#conventional amino acid & refuses to add H-atoms to it.

import os
import re
import chimera
from chimera import runCommand
path_unhadded_files = "D:\\unhadded_files\\"
path_hadded_files = "D:\\hadded_files\\"

if not os.path.isdir(path_hadded_files):
    os.makedirs(path_hadded_files)

for i in os.listdir(path_unhadded_files):
    if "+" in i:
        acceptor1 = i.split("+")            # this list will look something like "['1D5T_SD_MET_23_A', '_OE1_GLN_279_A.pdb']"
        acceptor2 = acceptor1[1].split("_") # this list will look something like "['', 'OE1', 'GLN', '279', 'A.pdb']"
        acceptor3 = ":" + acceptor2[3] + "." + acceptor2[4][0] + "@C" # it will look like ":279.A@C" & will be used to change the name of "C" atom of carbonyl carbon to "CA"
        acceptor4 = ":" + acceptor2[3] + "." + acceptor2[4][0]        # it will look like ":279.A" & will be used to change the name of acceptor residue to "AAA"
        acceptor5 = i[0:4] + acceptor1[1]
        acceptor6 = acceptor1[0] + ".pdb"

        #code for donor+acceptor file
        runCommand("open " + path_unhadded_files + i)
        if len(acceptor2[1]) == 1: #if the 2nd element of the list is equals 1, it means that the acceptor atom is a mainchain oxygen
            runCommand("setattr a name CA "  + acceptor3) # renaming the "C" atom of the acceptor carbonyl group to "CA" (reason explained in the beginning)
            runCommand("setattr r type AAA " + acceptor4) # renaming the acceptor aacid to "AAA" (reason explained in the beginning)
            runCommand("addh")
            new_name1 = i.split(".pdb")[0] + "_donaccep.pdb" # it will look like "1BKB_SG_CYS_30_A+_O_ASP_26_A_donaccep.pdb"
            runCommand("write 0 " + path_hadded_files + new_name1)
            runCommand("close 0")
        if len(acceptor2[1]) > 1: #if the 2nd element of the list is greater than 1, it means that the acceptor atom is a sidechain oxygen
            runCommand("addh")
            new_name1 = i.split(".pdb")[0] + "_donaccep.pdb"  # it will look like "1BKB_SG_CYS_30_A+_O_ASP_26_A_donaccep.pdb"
            runCommand("write 0 " + path_hadded_files + new_name1)
            runCommand("close 0")     
runCommand("stop")