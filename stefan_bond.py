# -*- coding: utf-8 -*-
"""
Created on Mon May 29 10:12:06 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import magnolia.needless_essential as ne
import sys
import matplotlib.pyplot as plt
import numpy as np

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_300_O2\Production\Sim-1'
filename    = 'bonds.reaxc'
bondfile_path    = directory+'\\'+filename

neighbours,atomtypes,bondorders,mtypes = bfp.get_neighbours(bondfile_path, bo=True, mtypes=True)
steps = list(neighbours.keys())
#%%--
timestep = 0.25
ps = bfp.step2picosecond(steps, timestep)
cutoff = 0.2
nobonds = np.array([0]*len(ps))
single  = np.array([0]*len(ps))
double  = np.array([0]*len(ps))
triple  = np.array([0]*len(ps))

for i,step_bo in enumerate(bondorders.items()):
    step,bo = step_bo
    for atom1, _ in bo.items():
        for atom2, bovalue in _.items():
            if bovalue<=cutoff:
                nobonds[i]+=1
            elif abs(bovalue-1)<=cutoff:
                single[i]+=1
            elif abs(bovalue-1.5)<=cutoff:
                double[i]+=1
            elif abs(bovalue-2.5)<=cutoff:
                triple[i]+=1

plt.plot(ps,nobonds/2,color='red',label='BO: 0 - {}'.format(cutoff))
plt.plot(ps,single/2,color='black',label='BO: {} - {}'.format(1-cutoff,1+cutoff))
plt.plot(ps,double/2,color='blue',label='BO: {} - {}'.format(1.5-cutoff,1.5+cutoff))
plt.plot(ps,triple/2,color='green',label='BO: {} - {}'.format(2.5-cutoff,2.5+cutoff))
plt.legend()