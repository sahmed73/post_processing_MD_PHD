# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 28 04:43:40 2024
"""

import pandas as pd
import dft_helper as dft
from rdkit.Chem import GetPeriodicTable

## user input
indirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\PAOr\pre'
outdirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\PAOr\post'

optfile=indirr+r'\gaussian_opt.out'
df=dft.extract_optimized_coordinates(optfile)

positions=df[['X','Y','Z']]
atomic_numbers=df['Atomic Number']

pt = GetPeriodicTable()
atomic_symbols = [pt.GetElementSymbol(num) for num in atomic_numbers]


write_infile_path=outdirr+r'\gaussian_opt.in'
dft.generate_gaussian_opt_input(write_infile_path, atomic_symbols, positions)
# dft.generate_xyzfile(outdirr+r'\opt_coords.xyz', atomic_symbols, positions)