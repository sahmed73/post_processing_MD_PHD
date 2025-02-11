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
from matplotlib.lines import Line2D

directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250\Sim-1"

filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bonddata  = bfp.parsebondfile(bondfile, charge=True)
#%%
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
charges    = bonddata['charge']
asyms      = ['H', 'C', 'O']

## filtering the carbons
def get_atomtypes(xxx):
    nd_atypes = []
    for x in xxx:
        nd_atypes.append(atypes[x])
    return np.array(nd_atypes)


seek_atomtypes = 3
for step, atom_charge in charges.items():
    temp_atom   = np.array(list(atom_charge.keys()))
    temp_charge = np.array(list(atom_charge.values()))
    mask   = (get_atomtypes(temp_atom) == seek_atomtypes) & (temp_atom <= 3200)
    atom   = temp_atom[mask]
    charge = temp_charge[mask]
    break # only one step

atom_ids = atom
charges = charge

# Reshape charges for KMeans
charges_reshaped = charges.reshape(-1, 1)

# Apply KMeans clustering
n_clusters = 2
kmeans = KMeans(n_clusters=n_clusters,random_state=42,n_init=20)
kmeans.fit(charges_reshaped)

# Get cluster assignments for each data point
labels = np.array(kmeans.labels_)

    
label_based_atoms = {}
for lab, atom in zip(labels,atom_ids):
    if lab not in label_based_atoms.keys():
        label_based_atoms[lab]=[atom]
    else:
        label_based_atoms[lab].append(atom)
        
label_based_charges = {}
for lab, charge in zip(labels,charge):
    if lab not in label_based_charges.keys():
        label_based_charges[lab]=[charge]
    else:
        label_based_charges[lab].append(charge)

plot_labels = []
for i in range(n_clusters):
    length = len(label_based_atoms[i])
    print("Label={}, Length={}, Mean charge={:.2f}".format(i,length,np.mean(label_based_charges[i])))
    

# Visualize the clusters
# Create scatter plots for each type
for i in range(n_clusters):
    plt.scatter(atom_ids[labels == i], charges[labels == i],s=7)
    
  
fsize = 20
plt.xlabel('{} atoms'.format(asyms[seek_atomtypes-1]),fontsize=fsize)
# plt.legend(ncol=3)
plt.xticks([])
plt.tick_params(axis='y', labelsize=fsize-5)
# plt.ylim(-0.6,0.6)
plt.ylabel('Charge',fontsize=fsize)
plt.title(directory[directory.find('ABCDE'):]+'\n'+'KMeans Clustering based on Charges\n', fontsize=8)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\Charge_Clusters_{}'.format(random.randint(0,999999999999999)), dpi=500,bbox_inches='tight')
plt.show()
    
