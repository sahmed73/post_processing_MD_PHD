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

baseoil = "PAO4"
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bonddata  = bfp.parsebondfile(bondfile, charge=True, abo=True)
#%%
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
charges    = bonddata['charge']
# abo        = bonddata['abo']

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

# Get cluster assignments for each data point
labels = np.array(kmeans.labels_)

# Analyze the clusters
for i in range(n_clusters):
    mean_charge = np.mean(charges[labels == i])
    print(f"Mean charge for cluster {i}: {mean_charge:.2f}")
    
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
colors = ['r','b','g','k','gold']
plt.style.use('default')
for i in range(n_clusters):
    plt.scatter(atom_ids[labels == i], charges[labels == i],s=15,
                       color = colors[i], alpha=0.3, linewidths=1)
    # plt.scatter(atom_ids[labels == i], charges[labels == i],s=15,
    #             facecolor=None,alpha=0.3,
    #             edgecolors='black', linewidths=1)
      

plt.xlabel('Carbon atoms')
# plt.legend(ncol=3)
plt.xticks([])
# plt.ylim(-0.36,0.01)
plt.ylabel('Charge')
plt.title(directory[directory.find('ABCDE'):]+'\n'+'KMeans Clustering based on Charges\n', fontsize=8)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\Charge_Clusters_{}'.format(random.randint(0,999999999999999)), dpi=500)
plt.show()
    
