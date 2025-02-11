# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Apr  4 13:43:14 2024
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO_Radical\20_PAO_Radical_20_A_Soria\Production\800K\Sim-1'
filename = '\\bonds.reaxc'
bondfile = dirr+filename
bonddata = bfp.parsebondfile(bondfile,cutoff=0.3,mtypes=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
steps     = np.array(list(neighbors.keys()))
ps = bfp.step2ps(steps,timestep=0.25)

rad  = np.array([92+x*92 for x in range(20)])
H_OH = np.array([1865+x*52 for x in range(20)])

Habs_counts = []
for step, neigh in neighbors.items():
    count = 0
    for parent, children in neigh.items():
        if parent in rad:
            for child in children:
                if atypes[child]==1 and child not in H_OH:
                    if 1<=mtypes[child]<=20:
                        count+=1
    Habs_counts.append(count)
plt.plot(ps,Habs_counts)
#%%
plt.style.use('classic')
plt_width = 7
aspect_ratio = 1.333333
plt.figure(figsize=[plt_width,plt_width/aspect_ratio])
plt.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.labelpad'] = 15
plt.rc('axes', grid=True)
plt.rc('font', size=10) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=20) # Axes label size
plt.rc('xtick', labelsize=20) # X-axis tick label size
plt.rc('ytick', labelsize=20) # Y-axis tick label size
plt.rc('legend', fontsize=20) # Legend fontsize

H_Sourcing_from = ['PAO', 'OH in AO', 'Other H in AO','NONE']
H_sourcing_count = [1,15,3,1]
plt.bar(H_Sourcing_from,H_sourcing_count)
plt.xticks(rotation=45)
plt.ylabel('H Counts')
plt.savefig('bar.png', dpi=500, transparent=True, bbox_inches='tight')