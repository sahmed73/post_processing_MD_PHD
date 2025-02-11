# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Oct 25 16:13:38 2023
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
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1'.format(baseoil,baseoil)
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bonddata  = bfp.parsebondfile(bondfile, bo=True)
#%%


neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
bondorders = bonddata['bondorders']
asyms      = ['H','C','O']

bo_list = []
for step, b_order in bondorders.items():
    for atom1,_ in b_order.items():
        for atom2, bo in _.items():
            bond = sorted((atom1,atom2))
            if atom1>atom2:
                bo_list.append(bo)
    break
                
bo_array = np.array(bo_list)


charges = bo_array
atom_ids = np.array(range(len(charges)))

# Reshape charges for KMeans
charges_reshaped = charges.reshape(-1, 1)

# Apply KMeans clustering
n_clusters = 2
kmeans = KMeans(n_clusters=n_clusters,n_init=10)
kmeans.fit(charges_reshaped)

# Get cluster assignments for each data point
labels = np.array(kmeans.labels_)
print(np.count_nonzero(labels==1))
# Analyze the clusters
colors = ['r','b','g','grey','orange','black']
for i in range(n_clusters):
    mean_charge = np.mean(charges[labels == i])
    min_charge  = np.min(charges[labels==i])
    max_charge  = np.max(charges[labels==i])
    print("Mean bond order for cluster {}: {:.2f}, min={}, max={}".format(i,mean_charge,min_charge,max_charge))
    plt.fill_between(range(np.max(atom_ids)),min_charge,max_charge,color=colors[i],alpha=0.15)
    
label_based_atoms = {}
for lab, atom in zip(labels,atom_ids):
    if lab not in label_based_atoms.keys():
        label_based_atoms[lab]=[atom]
    else:
        label_based_atoms[lab].append(atom)
        
    

# Visualize the clusters
# Create scatter plots for each type
for i in range(n_clusters):
    plt.scatter(atom_ids[labels == i], charges[labels==i],color=colors[i],s=7)
  

plt.xlabel('Bonds')
# plt.legend(ncol=3)
plt.xticks([])
plt.ylabel('Bond orders')
plt.title(directory[directory.find('Base'):]+'\n\n'+'KMeans Clustering based on bond orders\n', fontsize=8)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\BO_Clusters_{}'.format(random.randint(0,999999999999999)), dpi=500,bbox_inches='tight')
plt.show()
    
