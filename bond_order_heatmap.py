# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:14:53 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import pickle
import os
import time
import random
import seaborn as sns
import numpy as np
import pandas as pd

start_time = time.time()

timestep = 0.25
print('Assigned Timestep: ',timestep)

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_300_O2\Production\Sim-1'

filename    = 'bonds.reaxc'
bondfile    = directory+'\\'+filename

# check if the bond pickle file exists
picklefile     = directory+ '\\' + 'pickle_bond.pkl'
if os.path.exists(picklefile):
    print('Bondfile is loading from Pickle...')
    load_timestep, output_get_neighbours = pickle.load(open(picklefile, 'rb'))
    if load_timestep!=timestep:
        sys.exit('Error: Assigned timestep does not match with loaded the timestep\nAssigned Timestep: {} and Loaded Timestep: {}'.format(timestep,load_timestep))
else:
    print('Calling get_neighbours function....')
    output_get_neighbours = bfp.get_neighbours(bondfile,atypes=True,bo=True)
    with open(picklefile,'wb') as pf:
        pickle.dump((timestep,output_get_neighbours), pf)

neighbours,atomtypes,bondorders = output_get_neighbours
#%%

C_O_bond_3 = [(733, 734), (858, 859), (608, 609), (233, 234), (108, 109), (408, 409), (1108, 1109), (83, 84), (583, 584), (308, 309), (908, 909), (58, 59), (1233, 1234), (1183, 1184), (958, 959), (333, 334), (758, 759), (158, 159), (683, 684), (783, 784), (508, 509), (8, 9), (483, 484), (1158, 1159), (533, 534), (883, 884), (833, 834), (1208, 1209), (383, 384), (1033, 1034), (358, 359), (258, 259), (133, 134), (808, 809), (458, 459), (208, 209), (1058, 1059), (433, 434), (633, 634), (33, 34), (283, 284), (708, 709), (1008, 1009), (983, 984), (1083, 1084), (183, 184), (1133, 1134), (558, 559), (933, 934), (658, 659)]

bo_evolution = {}
steps = list(neighbours.keys())
ps = bfp.step2picosecond(steps, timestep)

for index, CO_ID in enumerate(C_O_bond_3):
    Cid,Oid = CO_ID
    bo_evolution[index+1]={}
    for step,p in zip(steps,ps):
        bo = bondorders[step][Cid].get(Oid)
        if bo==None: bo=0
        bo_evolution[index+1][p]=bo
        
        
df = pd.DataFrame(bo_evolution)
df = df.transpose()
df = df.loc[:,:1000.00]
df = df.sort_values(1000.00)
df.index = range(1,51)

#print(df)

_, ax1 = plt.subplots(figsize=(12,12))
cbar_kws = { 'ticks' : [0, 2], 'label': 'Bond order'}
hm = sns.heatmap(df,cmap='jet',cbar_kws={ 'label': 'Bond order'},xticklabels=1000,ax=ax1,vmin=0.0,vmax=3.0)
hm.set_xlabel('Time (ps)',fontsize=15)
hm.set_ylabel('Molecule',fontsize=15)
plt.yticks(rotation=0)
hm.set_title('Evolution of Bond Order for Type-3 C-O Bonds in all 50 AO1 Molecules',fontsize=15)


#colorbar settings
cbar = hm.collections[0].colorbar
#cbar.set_ticklabels([0.0,1.0,2.0,3.0])

fig = hm.get_figure()
savedir = 'python_outputs\\heatmaps\\heatmap'#+str(random.randint(100000, 999999))
#fig.savefig(savedir, dpi=400, bbox_inches='tight')







