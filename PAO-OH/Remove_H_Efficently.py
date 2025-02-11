# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Aug 29 00:38:08 2024


Calculate distance from all H to O in a single molecule
delete H whcih is close to the O
"""

import pandas as pd
import numpy as np


dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\1_ARLOI_24-08-28--00-09-01\15_PAO-OH_20_A0001\Equilibrate\Sim-1"

datafile = dirr+'\\equilibrated.npt.data'
outfile = dirr+'\\updated_equilibrated.npt.data'

natoms  = None
natypes = None

header     = ''
atoms      = None
velocities = None

start = None
end   = None

with open(datafile, 'r') as file:
    for i, line in enumerate(file):
        if start is None:
            header+=line
        if line.endswith('atoms\n'):
            natoms  = int(line.split()[0])
        if line.endswith('types\n'):
            natypes = int(line.split()[0])
            
        if line.startswith('Atoms'):
            start=i+2
        if line.startswith('Velocities'):
            end=i-2
    

atoms = pd.read_csv(datafile, delim_whitespace=True, skiprows=start,
                    nrows=end-start+1, header=None)

velocities = pd.read_csv(datafile, delim_whitespace=True,
                         skiprows=end+3, header=None)

#


