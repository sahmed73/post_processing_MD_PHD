# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Oct  1 03:39:14 2024
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


onsets = []
AOs    = ['A0001','A0002']#,'A0003', 'A0004','A0005']
AO_formula = dict(zip(AOs, ["H30C19O3", "H34C26O4", "H24C15O", "H24C17O2", "H44C29O2"]))

for i, AO in enumerate(AOs):
    dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\{}\Production\Sim-1".format(AO)
    bondfile = dirr+'\\bonds.reaxc'
    atomsymbols = ['H','C','O']   
    
    df = bfp.get_species_count(bondfile,atomsymbols,cutoff=0.3, timestep=0.25)
    
    x = 600+df['Time']*3.2941176470588234
    y = df[AO_formula[AO]]
    plt.scatter(x,y, alpha=0.05, label=AO)
    
    # onset, y_fit = bfp.get_onset(x, y, imc=15, ig=[1100,107,0,50])
    # onsets.append(onset)
    # plt.plot(x,y_fit, label=AO)
    plt.ylim(0, 15.9)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Antioxidant Count')
    plt.legend(loc="lower left")
#%%
    
#%%
fig, ax = plt.subplots(dpi=300)
ax.bar(AOs, onsets)
ax.set_xlabel("Antioxidant")
ax.set_ylabel("Onset of degradation temp (K)")

