# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 06:33:14 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\Production\1350\Sim-1'
bondfilepath = directory+'\\bonds.reaxc'

bondinfo = bfp.parsebondfile(bondfilepath)
#%%----
neighbours = bondinfo['neighbours']
atypes     = bondinfo['atypes']
atomsymbols= ['H','C','O']

swspecies  = bfp.get_SpeciesCountAtEveryTimestep(neighbours, atypes, atomsymbols)
#%%
def count_carbon(formula):
    i = formula.find('C')
    count = ''
    if i !=-1:
        for j in range(i+1,len(formula)):
            if formula[j].isdigit():
                count+=formula[j]
            if not formula[j].isdigit():
                break
        if not count: count+='1'
    if count:
        count = int(count)
    else:
        count = 0
    return count

n = 0
result = []
limit = 4500000
temp  = 1200
skipts= (1200-300)/4
print(limit*0.25/1000-skipts) #after this time count starts
for key in swspecies:
    carbon = [0]*16
    if key<limit: continue
    n+=1
    species = swspecies[key]
    for k in species:
        if k in ['O2','H62C30']: continue
        c = count_carbon(k)
        if c>15:
            carbon[-1]+=species[k]
        else:
            carbon[c]+=species[k]
    result.append(np.array(carbon))
carbon = sum(result)/n
std    = (sum([(x-carbon)**2 for x in result])/n)**0.5
x = [str(s) for s in range(len(carbon)-1)]+['>15']
plt.bar(x,carbon,yerr=std,color='maroon',capsize=5,alpha=0.9)
plt.xticks(fontsize=8)
plt.xlabel('Number of carbons')
plt.ylabel('Number of species')
plt.ylim(0)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\popo.png',dpi=400)
