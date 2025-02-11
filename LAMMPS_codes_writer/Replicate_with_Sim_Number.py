# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Apr  4 21:20:20 2024
"""

import os

base_dir = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\Pure_50B\Production\1050"

files = ['\Sim-1\input.in', '\Sim-1\submit.merced.sh', '\Sim-1\submit.sh']

number_of_simulation = 3

# Function to read file content
def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.read()

# Function to modify file content
def modify_content(file_content, sim_number):
    # Modify simulation number in input.in
    if 'variable\t\t\tsim_number equal' in file_content:
        return file_content.replace('variable\t\t\tsim_number equal 1', f'variable\t\t\tsim_number equal {sim_number}')
    # Modify job-name in submit files
    elif '#SBATCH --job-name=' in file_content:
        parts = file_content.split('\n')
        for i, line in enumerate(parts):
            if line.startswith('#SBATCH --job-name='):
                job_name = line.split('=')[1]
                new_job_name = job_name.replace('S1', f'S{sim_number}')
                parts[i] = '#SBATCH --job-name=' + new_job_name
                break
        return '\n'.join(parts)
    return file_content

# Main loop to create modified copies
for sim_number in range(2, number_of_simulation+1):
    sim_dir = os.path.join(base_dir, f'Sim-{sim_number}')
    os.makedirs(sim_dir, exist_ok=True)
    
    for file_path in files:
        file_path = base_dir + file_path
        content = read_file(file_path)
        modified_content = modify_content(content, sim_number)
        new_file_path = os.path.join(sim_dir, os.path.basename(file_path))
        
        with open(new_file_path, 'w') as new_file:
            new_file.write(modified_content)

