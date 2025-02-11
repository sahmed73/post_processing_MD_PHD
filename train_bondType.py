# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Nov 12 19:26:47 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import os
import sys
from sklearn.cluster import KMeans
import numpy as np

bond_orders = []
baseoil = "PAO4"
for temp in ['1200','1600']:
    for sim in ['Sim-1','Sim-2','Sim-3']:
        directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\{}\{}'.format(baseoil,baseoil,temp,sim)
        filename  = '\\bonds.reaxc'
        bondfile  = directory+filename
        print(temp,sim)
        bonddata = bfp.parsebondfile(bondfile,bo=True)
        neighbours = bonddata['neighbours']
        atypes     = bonddata['atypes']
        bondorder  = bonddata['bondorders']
        for step,innerDict_1 in bondorder.items():
            for _, innerDict in innerDict_1.items():
                for _, bo in innerDict.items():
                    if bo not in bond_orders:
                        bond_orders.append(bo)
#%%
bo_array = np.array(bond_orders).reshape(-1, 1)
# Apply KMeans clustering
n_clusters=3
kmeans = KMeans(n_clusters=n_clusters,n_init=10)
kmeans.fit(bo_array)

# Re-labeling the cluster. min_bo=Label-1, mid_bo=Label-2, max_bo=Label-3
# One problem. Only work if single-double-triple, this order maintained
# For example: single-triple = won't work; double-triple = won't work
# Only single = will work, single-double = will work
labels = kmeans.labels_
centroids = kmeans.cluster_centers_

# Pair up cluster labels with their centroids and sort
sorted_clusters = sorted(enumerate(centroids), key=lambda x: x[1])

# Create a mapping from old to new labels
label_mapping = {old_label: new_label+1 for new_label, (old_label, _) in enumerate(sorted_clusters)}

# Apply the new labels to your data
new_labels = [label_mapping[label] for label in labels]

# Visualize the clusters
# Create scatter plots for each type
colors = ['r','b','g','grey','orange','k','cyan']
# plt.figure(figsize=[20,15])
for i in range(n_clusters):
    plt.scatter(bo_array[labels==i],np.arange(bo_array.size)[labels == i],s=1,
                color=colors[i])
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures\bond_cluster',dpi=300,bbox_inches='tight')
plt.show()