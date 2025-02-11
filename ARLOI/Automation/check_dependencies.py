# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov  7 10:42:46 2024
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

def check_dependencies(dirr):
    # this will return if the dependencies exist
    input_file_path = os.path.join(dirr, 'input.in')
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
                        file_location = os.path.abspath(os.path.join(dirr, file_location))
                        print(f"Dependency: \"{file_location}\"")
                        
                    if os.path.exists(file_location):
                        return True
                    else:
                        return False


dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\test_delete_afterwards\A0005\Equilibration\Density_Eq\Test-1"
check_dependencies(dirr)