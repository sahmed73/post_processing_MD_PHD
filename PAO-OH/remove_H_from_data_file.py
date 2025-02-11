# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Jun 19 01:06:14 2024
"""
import pandas as pd
import re
import sys

dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0003\Equilibrate\Sim-1"
datafile = dirr+'\\equilibrated.npt.data'
outfile = dirr+'\\updated_equilibrated.npt.data'


#get the data in pandas dataframe

natoms  = None
natypes = None

header     = ''
atoms      = None
velocities = None

## get the start and end pointer of the atoms and velocities
start = None
end   = None

with open(datafile, 'r') as file:
    for i, line in enumerate(file):
        ## line.strip() wont work since its delete '/n'
        if start is None:
            ## save the headers until the word 'Atoms' is not found
            ## to understand please see the datafile format
            header+=line
        
        ## get the number of atoms
        if line.endswith('atoms\n'):
            natoms  = int(line.split()[0])
        
        ## get the number of types
        if line.endswith('types\n'):
            natypes = int(line.split()[0])
        
        ## get the line number where 'Atoms' located 
        if line.startswith('Atoms'):
            start=i+2 # 1-based line number, a blank line after 'Atoms'
        
        ## get the line number where all the atoms positon has finished
        if line.startswith('Velocities'):
            end=i-2 
    

atoms = pd.read_csv(datafile, delim_whitespace=True, skiprows=start,
                    nrows=end-start+1, header=None)

velocities = pd.read_csv(datafile, delim_whitespace=True,
                         skiprows=end+3, header=None)

##
# Hard Coding
N = 20  # How many leading molecules need H removed
r = 53  # remainder of H particle id upon division by 93 (PAO_radical)

# take evrything except the H atom that we dont want
updated_atoms      = atoms[~((atoms[0] % 93 == r) & (atoms[0] <= 93*N))].sort_values(by=0)
updated_velocities = velocities[~((velocities[0] % 93 == r) & (velocities[0] <= 93*N))].sort_values(by=0)

# updated number of atoms
updated_natoms = updated_atoms.shape[0]

updated_atoms[0]      = range(1, updated_natoms+1)
updated_velocities[0] = range(1, updated_natoms+1)

atom_count_regex = r'\d+ atoms'
atom_types_regex = r'\d+ atom types'

# Replace the matched values with the new values
updated_header = re.sub(atom_count_regex, f'{updated_natoms} atoms', header)
updated_header = re.sub(atom_types_regex, f'{natypes} atom types', updated_header)

with open(outfile, 'w') as f:
    f.write(updated_header)
    f.write('\n')

    # Write the updated atoms
    updated_atoms.to_string(f, index=False, header=False)
    f.write('\n\n')

    # Write the velocities section header
    f.write('Velocities\n\n')
    updated_velocities.to_string(f, index=False, header=False)
    f.write('\n')