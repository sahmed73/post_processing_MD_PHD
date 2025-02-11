# -*- coding: utf-8 -*-
"""
Created on Wed May 10 01:54:29 2023

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
#%%----

C_O_bond_4 = [(613, 612), (88, 87), (863, 862), (413, 412), (1113, 1112), (238, 237), (113, 112), (1238, 1237), (1188, 1187), (313, 312), (588, 587), (263, 262), (913, 912), (63, 62), (338, 337), (688, 687), (1038, 1037), (963, 962), (513, 512), (213, 212), (788, 787), (438, 437), (763, 762), (638, 637), (1088, 1087), (838, 837), (1163, 1162), (163, 162), (888, 887), (538, 537), (563, 562), (388, 387), (738, 737), (188, 187), (363, 362), (813, 812), (38, 37), (488, 487), (138, 137), (1013, 1012), (713, 712), (1138, 1137), (988, 987), (13, 12), (463, 462), (288, 287), (663, 662), (1063, 1062), (1213, 1212), (938, 937)]

steps = list(neighbours.keys())
y = [0]*len(steps)
for atomid_to_be_tracked, neighbouring_carbon in C_O_bond_4:  
    for index,step in enumerate(steps):
        bo = bondorders[step][atomid_to_be_tracked].get(neighbouring_carbon)
        if bo==None: bo=0
        y[index]+=bo
    

x = np.array(bfp.step2picosecond(steps,timestep))
y = np.array(y)/len(C_O_bond_4)

atomsymbols = 'HCO' 
plt.style.use('seaborn')
plt.scatter(x,y,s=100, alpha=0.6, edgecolor='black', linewidth=1)
plt.xlabel('Time (ps)')
plt.ylabel('Bond Order')
plt.title('Atom {} ({}) and {} ({}) bond order evolution'.format(atomid_to_be_tracked,atomsymbols[atomtypes[atomid_to_be_tracked]-1],neighbouring_carbon,atomsymbols[atomtypes[neighbouring_carbon]-1]))
savedir = 'python_outputs\\figures\\bond_order'
plt.savefig(savedir, dpi=300, bbox_inches='tight')