# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:36:18 2023

@author: arup2
"""

import magnolia.dumpfile_parser as dfp
import pandas as pd
import matplotlib.pyplot as plt


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_200_O2\Production\Sim-1'
filename = '\\oxidation.lammpstrj'
dumpfile = directory+filename

dumpdata = dfp.parsedumpfile(dumpfile)
dist     = dfp.distance_tracker(dumpdata,480,483)

plt.style.use('seaborn')
plt.scatter(dist.keys(),dist.values(),s=100, alpha=0.6, edgecolor='black', linewidth=1)