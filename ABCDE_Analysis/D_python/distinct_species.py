# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Mar 27 18:15:08 2024
"""
import magnolia.speciesfile_parser as sfp
import pandas as pd

dirrs = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\Sim-1']

filename  = '\\species.out'

for dirr in dirrs:
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile)
    distinct_count = species.index.size