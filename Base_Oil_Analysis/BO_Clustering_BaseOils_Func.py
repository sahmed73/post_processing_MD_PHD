# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Oct 25 20:48:25 2023
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import numpy as np
import random
from sklearn.cluster import KMeans
from PIL import Image, ImageSequence
import os


def bo_clustering(bondorders,n_clusters,atypes,frame=0,fill=False):
    bo_list = []
    steps = list(bondorders.keys())
    for atom1,bo_ in bondorders[steps[frame]].items():
        if atypes[atom1]!=2: continue
        for atom2, bo in bo_.items():
            if atypes[atom2]!=2: continue
            if atom1>atom2:
                bo_list.append(bo)
                    
    bo_array = np.array(bo_list)


    charges = bo_array
    atom_ids = np.array(range(len(charges)))

    # Reshape charges for KMeans
    charges_reshaped = charges.reshape(-1, 1)

    # Apply KMeans clustering
    kmeans = KMeans(n_clusters=n_clusters,n_init=10)
    kmeans.fit(charges_reshaped)

    # Get cluster assignments for each data point
    cluster_labels = np.array(kmeans.labels_)
    labels = cluster_labels
    
    # Custom Labels
    # labels_meanCharge = []
    # for i in range(n_clusters):
    #     labels_meanCharge.append(np.mean(charges[cluster_labels == i]))
        
    # mean_charge  =  np.array([])
    # for lab in cluster_labels:
    #     mean_charge = np.append(mean_charge, labels_meanCharge[lab])

    # deg2charge = np.array([1.,1.5]) # deg 1, 2, 3
    # labels = np.array([],dtype=int)
    # for each_charge in mean_charge:
    #     error    = np.abs(deg2charge-each_charge)
    #     minErrID = error.argmin()
    #     labels    = np.append(labels,minErrID)
    
    # Analyze the clusters
    colors = ['r','b','g','grey','orange','black']
    fig, ax = plt.subplots()
    for i in range(n_clusters):
        mean_charge = np.mean(charges[labels == i])
        min_charge  = np.min(charges[labels==i])
        max_charge  = np.max(charges[labels==i])
        print("Mean bond order for cluster {}: {:.2f}, min={}, max={}".format(i,mean_charge,min_charge,max_charge))
        if fill:
            ax.fill_between(range(np.max(atom_ids)),min_charge,max_charge,color=colors[i],alpha=0.15)        

    # Visualize the clusters
    # Create scatter plots for each type
    for i in range(n_clusters):
        ax.scatter(atom_ids[labels == i], charges[labels==i],color=colors[i],s=7)
           
    return labels,fig,ax
    
    

#%%
baseoil = "PAO4"
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1'.format(baseoil,baseoil)
filename  = '\\bonds.reaxc'
bondfilepath  = directory+filename

bonddata  = bfp.parsebondfile(bondfilepath, bo=True)
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
bondorders = bonddata['bondorders']
#%%
steps = np.array(list(bondorders.keys()))
dirr =  r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'

for frame in range(0,steps.size,10):
    labels, fig, ax = bo_clustering(bondorders, 3,atypes,fill=True)
    ax.set_xlabel('Bonds')
    ax.set_xticks([])
    ax.set_ylabel('Bond orders')
    # ax.set_ylim(0,3)
    # ax.text(0,2.85,"Double Bond: {}".format(list(labels).count(1)),color='b',fontsize=8)
    # ax.text(0,2.70,"Single Bond: {}".format(list(labels).count(0)),color='r',fontsize=8)
    # ax.text(0,2.55,"Total Bond: {}".format(list(labels).count(1)+list(labels).count(0)),
    #         color='k',fontsize=9)
    ax.set_title('frame: {}'.format(frame))
    fname = dirr+'\\BO_Clusters_{}.png'.format(frame)
    # fig.savefig(fname, dpi=500,bbox_inches='tight')
    plt.show()
    break