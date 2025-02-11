# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Sep  2 03:09:56 2024
"""

# H location

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0002\Production\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
speciesfile = directory+'species.out'
atomsymbols = ['H','C','O']
atominfo = bfp.parsebondfile(bondfile,bo=True)
#%%
neighbours = atominfo['neighbours']
atypes = atominfo['atypes']
bondorders = atominfo['bondorders']
atomsys = 'HCO'

species = sfp.get_species_count(speciesfile)
