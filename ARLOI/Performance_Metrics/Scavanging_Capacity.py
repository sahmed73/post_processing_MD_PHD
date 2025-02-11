# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri May 24 10:01:20 2024
"""

import magnolia.speciesfile_parser as sfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('default')
plt.rcParams['font.size'] = 18

molecules=["A","B",'BHT']
capacity = []
fig, ax = plt.subplots(dpi=350)
for molecule in molecules: 
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\20PAOtcr_15{}\Production\TempRamp'.format(molecule)
    filename = '\\species.out'
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile, timestep=0.25, Nfreq = 1)
    radicals = species['H61C30O']
    cap = (radicals.iloc[0]-radicals.iloc[-1])/15
    capacity.append(cap)
    
ax.bar(molecules,capacity,color='tab:green', edgecolor='k')
ax.set_ylabel("Scavenging Capacity")
    