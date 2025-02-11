# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 23:57:36 2023

@author: Shihab 
"""

import os
import shutil
import re
import io


def modify_file_content(filename, changes):
    with open(filename, 'r') as file:
        content = file.read()
        for old_text, new_text in changes.items():
            content = content.replace(old_text, new_text)
    with open(filename, 'w') as file:
        file.write(content)

def replace_line_with_string(file_path, search_string, replacement):
    """
    Replace lines in a file that contain a specific string.

    Parameters:
    - file_path (str): The path to the text file.
    - search_string (str): The string to search for in each line.
    - replacement (str): The string to replace the entire line with if search_string is found.
    """

    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify lines containing the search_string
    with open(file_path, 'w') as file:
        for line in lines:
            if search_string in line:
                file.write(replacement + "\n")
            else:
                file.write(line)

def convert_line_endings(filename):
    with io.open(filename, 'r', newline='\r\n', encoding='utf-8') as file:
        lines = file.readlines()
    
    with io.open(filename, 'w', newline='\n', encoding='utf-8') as file:
        file.writelines(lines)


        
def generate_sim(common_dir, simpath="\\Sim-1", Nsim=1):
    source_dir  = common_dir+simpath
    for i in range(Nsim-1):
        sim_number = i+2
        dest_dir    = common_dir+"\\Sim-{}".format(sim_number)
        # Ensure destination doesn't exist, then copy
        if os.path.exists(dest_dir):
            print('Directory already exist: ',dest_dir)
            break
        shutil.copytree(source_dir, dest_dir)
    
        # Now let's modify the files in the copied directory
        files_to_modify = {
            "input.in": {
                "sim_number equal 1": "sim_number equal {}".format(sim_number),
                # Add more replacements as needed
            },
            "submit.merced.sh": {
                "job-name=S1-1600_PAO": "job-name=S{}-1600_PAO".format(sim_number),
                # Add more replacements as needed
            },
            "submit.sh": {
                "job-name=S1-1600_PAO": "job-name=S{}-1600_PAO".format(sim_number),
                # Add more replacements as needed
            }
        }
        for filename, changes in files_to_modify.items():
            modify_file_content(os.path.join(dest_dir, filename), changes)
            
        
def generate_temp(common_dir, source_temp, templist):
 
    source_dir = os.path.join(common_dir, str(source_temp))
    for temp in templist:
        dest_dir    = os.path.join(common_dir, str(temp))
        
        # Ensure destination doesn't exist, then copy
        if os.path.exists(dest_dir):
            print('Directory already exist: ',dest_dir)
            continue
        shutil.copytree(source_dir, dest_dir)
        
        # Now let's modify the files in the copied directory
        for dirpath, dirnames, filenames in os.walk(dest_dir):
            for dirname in dirnames:
                simpath   = os.path.join(dirpath, dirname)
                
                
                ## input file
                filepath = os.path.join(simpath,'input.in')
                search = "T equal"
                replacement = "variable   			T equal {}".format(temp)
                replace_line_with_string(filepath, search, replacement)
                # converting \r\n to \n
                convert_line_endings(filepath)
                    
                ## submit.sh file
                filepath = os.path.join(simpath,'submit.sh')
                with open(filepath,'r') as file:
                    search = "#SBATCH --job-name"  
                    for line in file:
                        if search in line:
                            jobname = line.strip()
                            break
                    # replace numbers that are surrounded by '-' or '_'
                    pattern = r'(?<=[-_])\d+(?=[-_])'
                    replacement = re.sub(pattern, str(temp), jobname)
                    replace_line_with_string(filepath, search, replacement)
                # converting \r\n to \n
                convert_line_endings(filepath)
                
                ## submit.merced.sh file
                filepath = os.path.join(simpath,'submit.merced.sh')
                with open(filepath,'r') as file:
                    search = "#SBATCH --job-name"  
                    for line in file:
                        if search in line:
                            jobname = line.strip()
                            break
                    # replace numbers that are surrounded by '-' or '_'
                    pattern = r'(?<=[-_])\d+(?=[-_])'
                    replacement = re.sub(pattern, str(temp), jobname)
                    replace_line_with_string(filepath, search, replacement)
                # converting \r\n to \n
                convert_line_endings(filepath)
                    
                
                    
        
        

        
common_dir = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\25_Squalane_200_O2_Soria\Production"
# simpath = "\\Sim-1"

# generate_sim(common_dir,simpath,Nsim=3)
generate_temp(common_dir,1200,templist=[1500,1550,1600])



