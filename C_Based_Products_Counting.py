# -*- coding: utf-8 -*-
"""
Created on Tue May 16 13:20:52 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_300_O2\Production\Sim-1'

filename  = '\\bonds.reaxc'
bondfile  = directory+filename

neighbours,atomtypes = bfp.get_neighbours_v2(bondfile)

#%%

steps       = list(neighbours.keys())
ps          = bfp.step2picosecond(steps, 0.25)

count_C0_3  = []
count_C4_6  = []
count_C7_9  = []
count_C10   = []

carbon_type = 2 # C 

for step, neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    temp_count_1 = 0
    temp_count_2 = 0
    temp_count_3 = 0
    temp_count_4 = 0
    for molecule in molecules:
        mol_with_type = list(map(lambda x:atomtypes[x],molecule))
        cc = mol_with_type.count(carbon_type)
        if 0<=cc<=3:
            temp_count_1 += 1
        elif 4<=cc<=6:
            temp_count_2 +=1
        elif 7<=cc<=9:
            temp_count_3 +=1
        elif cc==10:
            temp_count_4 +=1
    count_C0_3.append(temp_count_1)
    count_C4_6.append(temp_count_2)
    count_C7_9.append(temp_count_3)
    count_C10.append(temp_count_4)
#%%
plt.plot(ps,count_C0_3,color='black')
plt.plot(ps,count_C4_6,color='blue')
plt.plot(ps,count_C7_9,color='red')
plt.plot(ps,count_C10,color='green')
plt.xlabel('Time (ps)')
plt.ylabel('Molecule Count')
plt.savefig('C_Bases_Products_4', dpi=400)