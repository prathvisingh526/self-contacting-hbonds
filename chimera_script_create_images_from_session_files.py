import os
from chimera import runCommand

path_session_files = "D:\\BBBB\\"
for file in os.listdir(path_session_files):
    if file[-3:] == ".py":
        filename_jpg = path_session_files + file.split(".")[0] + ".jpg" # ensure that "file" variable has no whitespaces, else the code will give error
        runCommand("open " + path_session_files + file)
        runCommand("copy file " + filename_jpg + " jpeg supersample 3")
        runCommand("close session")

for file in os.listdir(path_session_files):
    if file[-4:] == ".pyc":
        os.remove(path_session_files + file)
runCommand("stop") # exits chimera

