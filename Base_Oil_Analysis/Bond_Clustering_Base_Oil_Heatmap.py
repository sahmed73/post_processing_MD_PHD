# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 17:06:16 2023
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import numpy as np
import seaborn as sns
import random
from sklearn.cluster import KMeans

baseoil = "PAO4"
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1'.format(baseoil,baseoil)
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bonddata  = bfp.parsebondfile(bondfile, ALL=True)
#%%
neighbours = bonddata['neighbours']
mtypes     = bonddata['mtypes']
atypes     = bonddata['atypes']
charges    = bonddata['charge']
NLP        = bonddata['nlp']
bondorders = bonddata['bondorders']
firststepneigh = neighbours[list(neighbours.keys())[0]]

## filtering the data
def get_atomtypes(xxx):
    nd_atypes = []
    for x in xxx:
        nd_atypes.append(atypes[x])
    return np.array(nd_atypes)

for step, atom_charge in charges.items():
    temp_atom   = np.array(list(atom_charge.keys()))
    temp_charge = np.array(list(atom_charge.values()))
    mask   = get_atomtypes(temp_atom) == 2
    atom   = temp_atom[mask]
    charge = temp_charge[mask]
    break # only one step

atom_ids = atom
charges = charge

# Reshape charges for KMeans
charges_reshaped = charges.reshape(-1, 1)

# Apply KMeans clustering
n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(charges_reshaped)

# Custom label based on carbon degree
cluster_labels = kmeans.labels_
labels_meanCharge = []
for i in range(n_clusters):
    labels_meanCharge.append(np.mean(charges[cluster_labels == i]))

mean_charge  =  np.array([])
for lab in cluster_labels:
    mean_charge = np.append(mean_charge, labels_meanCharge[lab])

deg2charge = np.array([-0.30,-0.20,-0.05]) # deg 1, 2, 3
label = np.array([],dtype=int)
for each_charge in mean_charge:
    error    = np.abs(deg2charge-each_charge)
    minErrID = error.argmin()
    label    = np.append(label,minErrID)

atom2label = {}
for lab,atom in zip(label,atom_ids):
    atom2label[atom]=lab
####################-Bonds-##############################
bonds = {}
for parent, children in firststepneigh.items():
    if atypes[parent] != 2:
        continue
    parent_label = atom2label[parent]
    for child in children:
        if atypes[child] != 2:
            continue
        child_label = atom2label[child]
        bond_label  = "".join(sorted(str(parent_label+1)+str(child_label+1)))
        
        bond = tuple(sorted((parent,child)))
        if bond_label in bonds:
            bonds[bond_label].add(bond)
        else:
            bonds[bond_label] = set()
print(len(bonds.keys()))
for key, value in bonds.items():
    print(key,value)
#%%####
## User input
plot_keys     = sorted(bonds.keys())
skipts        = (1600-300)/4
timestep      = 0.25
steps         = np.array(list(bondorders.keys()))
cbar_min      = 0.0
cbar_max      = 2.0
fontsize      = 15 # All fonts and labels
ticksize      = 10 # Only x-ticks and cbar ticks

# each plot use figsize=(2.25,4)
n_plot = len(plot_keys)
fig, ax = plt.subplots(1,len(plot_keys)+1,figsize=(2.25*n_plot, 4),
                       gridspec_kw={'width_ratios': [*[1]*n_plot, 0.1]})

for i,key in enumerate(plot_keys):
    b = bonds[key]
    df = bfp.bondorder_evolution(bondorders, b, ts=timestep,skipts=skipts,plot='no')

    if i == n_plot - 1:
        cbar_ticks = range(int(cbar_min),int(cbar_max)+1)
        cbar_kws = {'ticks': cbar_ticks}
        hm = sns.heatmap(df, cmap='jet', xticklabels=1000, ax=ax[i], 
                         vmin=cbar_min, vmax=cbar_max, cbar=True,
                         cbar_ax=ax[-1],cbar_kws=cbar_kws)
    else:
        hm = sns.heatmap(df, cmap='jet', xticklabels=1000, ax=ax[i], 
                         vmin=cbar_min, vmax=cbar_max, cbar=False)
    
    xticklabels = [int(float(x.get_text())) for x in ax[i].get_xticklabels()]
    ax[i].set_xticklabels(xticklabels,rotation=90,fontsize=ticksize)
    ax[i].set_yticks([])
    ax[i].set_title('{}'.format(key))

# Common X Label
fig.text(0.5, -0.04, 'Time (ps)', ha='center', va='center',fontsize=fontsize)
## Accessing the color bar. we created colorbar on our own. So cbar = ax[-1]
ax[-1].set_ylabel("Bond order", fontsize=fontsize, labelpad=10)
ax[-1].tick_params(labelsize=ticksize)

title = directory[directory.find('Base'):]
fig.suptitle(title, y=1.04, fontsize=10)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution\bo_evolution_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
