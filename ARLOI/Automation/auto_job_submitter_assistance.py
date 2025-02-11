# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov  5 11:25:06 2024

I wanted to prioriize the simulations but failed right now.
I will do it later. But now it will run simulation simply one-by-one
"""

import os
import re

def extract_lammps_variables(input_file_path):   
    variables = {}
    # this styels are from lampps documantation
    styles = "(delete|atomfile|file|format|getenv|index|internal|loop|python|string|timer|uloop|universe|world|equal|vector|atom)"
    variable_pattern = re.compile(rf"^\s*variable\s+(\w+)\s+{styles}\s+(.*)", re.IGNORECASE)

    with open(input_file_path, 'r') as file:
        for line in file:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            match = variable_pattern.match(line)
            if match:
                var_name = match.group(1) 
                var_value = match.group(3).strip() 
                variables[var_name] = var_value
    
    return variables

def find_simulation_dirr(root_dir):
    # this will return the simulation directories and 
    # dependent file path
    dirr = {}
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if 'input.in' in filenames:
            input_file_path = os.path.join(dirpath, 'input.in')
            variables = extract_lammps_variables(input_file_path)
            
            with open(input_file_path, 'r') as file:
                for line in file:                    
                    if line.strip().startswith('read_restart') or line.strip().startswith('read_data'):
                        parts = line.split()
                        if len(parts) > 1:
                            file_location = parts[1].strip()
                            pattern = r"\$\{?(\w+)\}?"
                            while True:
                                match = re.search(pattern, file_location)
                                if match:
                                    var_name = match.group(1)
                                    var_value = variables.get(var_name, f"${{{var_name}}}")
                                    file_location = re.sub(pattern, var_value, file_location, count=1)
                                else:
                                    break
                                

                            if not os.path.isabs(file_location):
                                file_location = os.path.abspath(os.path.join(dirpath, file_location))
                                
                            if os.path.exists(file_location):
                                dirr[dirpath]=(file_location,'exist')
                            else:
                                dirr[dirpath]=(file_location,'not_exist')
                            break

    return dirr

def prioritize_simulation_simple():
    '''
        |---Equilibration
        |    	|---Density_Eq
        |    	|		|---Sim-1
        |   	|---Energy_Eq
        |   			|---Sim-1
        |   			|---Sim-2
        |    			|---Sim-3
        |    			|---Sim-4
        |    			|---Sim-5
        |    			|---Sim-6
        |    			|---Sim-7
        |    			|---Sim-8
        |    			|---Sim-9
        |    			|---Sim-10
        |---Production
        |	    |---300-500K_TRR=1Kpps
        |		    	|---Sim-1
        |			    |---Sim-2
        |			    |---Sim-3
        |			    |---Sim-4
        |			    |---Sim-5
        |			    |---Sim-6
        |			    |---Sim-7
        |			    |---Sim-8
        |			    |---Sim-9
        |			    |---Sim-10
'''
    pass
    


root_dir = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001"
directories = find_simulation_dirr(root_dir)

output_file = root_dir+"\\simulation_directories.txt"
with open(output_file, 'w') as file:
    for dirr in directories:
        file.write(f"{dirr}\n{directories[dirr]}\n\n")