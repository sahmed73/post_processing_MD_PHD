# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov  3 10:58:23 2023
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
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\100_Less_TempRamp'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bonddata  = bfp.parsebondfile(bondfile, charge=True)
#%%


neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
charges    = bonddata['charge']


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
colors = ['#f7a0da', '#8f8f8f', '#86c7b1'] # pink, grey, cyan
for i in range(n_clusters):
    mean_charge = np.mean(charges[labels == i])
    min_charge  = np.min(charges[labels==i])
    max_charge  = np.max(charges[labels==i])
    print("Mean bond order for cluster {}: {:.2f}, min={}, max={}".format(i,mean_charge,min_charge,max_charge))
    plt.fill_between(range(np.max(atom_ids)),min_charge,max_charge,color=colors[i],alpha=0.3)
    
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

cluster = {"PAO4": {100: "Primary",600: "Secondary",50: "Tertiary"},
           "Squalane": {200: "Primary",400: "Secondary",150: "Tertiary"}}

plot_labels = []
for i in range(n_clusters):
    length = len(label_based_atoms[i])
    print("Label={}, Mean charge={:.2f}, Length={}, Cluster={}".format(i,np.mean(label_based_charges[i]),length,cluster[baseoil][length]))
    print(label_based_atoms[i])
    print('------------------------')
    plot_labels.append(cluster[baseoil][length])
    

# Visualize the clusters
# Create scatter plots for each type

for i in range(n_clusters):
    plt.scatter(atom_ids[labels == i], charges[labels == i], color=colors[i],s=7)
    
# plt.scatter(atom_ids, charges, c=labels, cmap='rainbow',s=7,marker='o')
# Custom legend
# for i in range(3):
#     plt.scatter([], [],color=colors[i],s=55, label="{} C".format(plot_labels[i]))
  

plt.xlabel('Carbon atoms')
# plt.legend(ncol=3)
plt.xticks([])
plt.ylim(-0.39,-0.001)
plt.ylabel('Charge')
plt.title(directory[directory.find('Base'):]+'\n\n'+'KMeans Clustering based on Charges\n', fontsize=8)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\Charge_Clusters_{}'.format(random.randint(0,999999999999999)), dpi=500,bbox_inches='tight')
plt.show()
    
