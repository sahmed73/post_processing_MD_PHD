# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 28 16:10:13 2024
"""

import pandas as pd
import dft_helper as dft
from rdkit.Chem import GetPeriodicTable
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

## user input
dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\A0005'

optfile=dirr+r'\opt.out'
df=dft.extract_optimized_coordinates(optfile)

positions=df[['X','Y','Z']]
atomic_numbers=df['Atomic Number']

pt = GetPeriodicTable()
atomic_symbols = [pt.GetElementSymbol(num) for num in atomic_numbers]

dft.write_xyzfile(optfile+'.xyz', atomic_symbols, positions)