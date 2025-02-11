# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri May 24 02:08:26 2024
"""
import sys
import re
import magnolia.speciesfile_parser as sfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('default')
plt.rcParams['font.size'] = 18

AOs = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']
fig, ax = plt.subplots(dpi=350)
for AO in AOs: 
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_{}\Production\Sim-1'.format(AO)
    if AO=='A0001':
        dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\1_ARLOI_24-08-28--00-09-01\15_PAO-OH_20_A0001\Production\Sim-2_with_H_removed'
    filename = '\\species.out'
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile, timestep=0.25, Nfreq = 1)
    
    radicals = species['H61C30O']
    time = species['Time']/1000 # ns
    
    ## data stripping
    N_values = [4,5,6,7,8]
    m_values = []
    for N_scav in N_values:
        N = radicals.max()-N_scav
        print(N,N_scav)
        idx=None if N is None else (radicals[::-1]==N).idxmax()+1
        radicals=radicals.loc[:idx]
        time=time.loc[:idx]
        log_radicals = np.log(radicals)
        
        ## fitting
        m,c  = np.polyfit(time, log_radicals, 1)
        m_values.append(-m)
    
    ax.plot(N_values,m_values,'o-', label=AO)
    ax.set_xlabel("Number of scavanged radicals")
    ax.set_ylabel("$k_{scav} \ (ns^{-1})$")
    # ax.set_title(f"Antioxidant {molecule}", fontsize=16)
    ax.legend(fontsize=12)
    ax.set_ylim(0.75, 13.9)

# ax.text(0.45,0.9,f"$k_{{scav}}={-m:.4} \ ns^{{-1}}$",transform=ax.transAxes, fontsize=16)