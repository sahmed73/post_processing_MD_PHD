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

baseoil = "Squalane"
sim_dir = ['Sim-1','Sim-2','Sim-3']
bonds = {}
for sim in sim_dir:
    directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}'.format(baseoil,baseoil,sim)
    filename  = '\\bonds.reaxc'
    bondfile  = directory+filename
    bonddata  = bfp.parsebondfile(bondfile, bo=True, charge=True)
    neighbours = bonddata['neighbours']
    atypes     = bonddata['atypes']
    charges    = bonddata['charge']
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
    local_bonds = {}
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
            if bond_label in local_bonds:
                local_bonds[bond_label].add(bond)
            else:
                local_bonds[bond_label] = set()
    
    bonds |= local_bonds
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
fig, ax = plt.subplots(1,len(plot_keys),figsize=(4.25*n_plot, 4))

stat_reactivity = {}
for i,key in enumerate(plot_keys):
    b = bonds[key]
    df = bfp.bondorder_evolution(bondorders, b, ts=timestep,skipts=skipts,plot='no')
    df_bool = (df<1.1) & (df>0.45)
    df_count= df_bool.sum(axis=0)
    # df_norm = df_count
    df_norm = (df_count-df_count.min())/(df_count.max()-df_count.min())
    
    # Fitting a linear function
    x = df_norm.index.values
    y = df_norm.values
    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit)
    
    stat_reactivity[key] = abs(10000*fit[0]) #taking the slope
    print(fit[0]*1000)
    
    # Plotting the data and the fit
    ax[i].scatter(x, y, s=5,label='simulation')
    ax[i].plot(x, fit_fn(x), color='red',label='fit')
    ax[i].set_title('{}'.format(key))
    ax[i].legend()

# Common X Label
fig.text(0.5, -0.04, 'Time (ps)', ha='center', va='center',fontsize=fontsize)

title = directory[directory.find('Base'):]
fig.suptitle(title, y=1.04, fontsize=10)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\bondplot'

plt.savefig(savedir+'\\bo_count_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
plt.show()
plt.bar(stat_reactivity.keys(),stat_reactivity.values())
plt.xlabel('Types of bonds',fontsize=fontsize)
plt.ylabel('Statistical bond transition rate ($x 10^4$)',fontsize=fontsize)
plt.savefig(savedir+'\\bond_transition_rate_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
