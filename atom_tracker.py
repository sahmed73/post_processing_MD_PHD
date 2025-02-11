# -*- coding: utf-8 -*-
"""
Created on Tue May  9 18:37:27 2023

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

atomsymbols = 'HCO'
atomid_to_be_tracked = 487
neighbouring_carbon  = 488
initial_neighbours    = neighbours[0][atomid_to_be_tracked]
print(initial_neighbours)
print(list(map(lambda x: atomtypes[x],initial_neighbours)))

steps = list(neighbours.keys())
y = []
x = []
flag = True
for step in steps:
    bo = bondorders[step][atomid_to_be_tracked].get(neighbouring_carbon)
    #print(bo,end='\t')
    if bo==None:
        y.append(0)
        if flag: 
            print('Bond Breaks at timestep: ',step)
            flag = False
    else:
        y.append(bo)
    x.append(step)

x = bfp.step2picosecond(x,timestep)



plt.style.use('seaborn')
plt.scatter(x,y,s=100, alpha=0.6, edgecolor='black', linewidth=1)
plt.xlabel('Time (ps)')
plt.ylabel('Bond Order')
plt.title('Atom {} ({}) and {} ({}) bond order evolution'.format(atomid_to_be_tracked,atomsymbols[atomtypes[atomid_to_be_tracked]-1],neighbouring_carbon,atomsymbols[atomtypes[neighbouring_carbon]-1]))
savedir = 'python_outputs\\figures\\bond_order'
plt.savefig(savedir, dpi=300, bbox_inches='tight')











