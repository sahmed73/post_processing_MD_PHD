# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Mar  7 00:48:58 2024
"""

import pandas as pd

#######-User Input-################################################

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\OPLSAA\PAO+Molecule-A'
filename  = r'PAO4.pdb'
columns = [2,-3,-2,-1] # columns of atom_name, x, y, z

###################################################################


pdbfile   = directory+'\\'+filename
xyzfile   = pdbfile[:pdbfile.rfind('.')]+'.xyz'

pdbdata = pd.read_csv(pdbfile, header=None, delim_whitespace=True, skiprows=1)

## xyz format
# taking only the rows that have 'ATOM' as their first entry
xyzdata = pdbdata[pdbdata[0]=='ATOM'].iloc[:,columns]
# taking the first character only for the atom symbol ## BE CAREFUL!!
xyzdata.iloc[:,0] = xyzdata.iloc[:,0].apply(lambda x:x[0])

# write in dot-aligned format. to_string preserve the dot-aligned format
string = xyzdata.to_string(index=False,header=False)
catch  = filename[:filename.rfind('.')] # optional. But must for fftool
second_line = f'{catch} {catch}.ff\n'

with open(xyzfile, "w") as file:
    # first line = number of atoms
    file.write(f'{xyzdata.shape[0]}\n')
    # second line = ff file path for fftool. Otherwise optional
    file.write(second_line)
    # writing the remainings
    file.write(string)

