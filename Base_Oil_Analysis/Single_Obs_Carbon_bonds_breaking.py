# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Dec  4 14:07:50 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import pandas as pd
import seaborn as sns
import numpy as np

baseoil = "PAO4"
sim_dirs = ['Sim-1','Sim-2','Sim-3']
counts = {x:0 for x in range(1,30)}
for sim in sim_dirs:
    directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
    bondfilepath = directory+'\\bonds.reaxc'
    cutoff = 0.45
    bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
    
    neighbors = bonddata['neighbours']
    atypes    = bonddata['atypes']
    mtypes    = bonddata['mtypes']
    bondorders= bonddata['bondorders']
    atomsymbols = 'HCO'
    
    baseoil_molecules = []  # no hydrogen
    for step, neigh in neighbors.items():
        main_graph = nx.Graph(neigh)
        molecules = bfp.get_molecules(neigh)
        for molecule in molecules:
            if len(molecule)==92:
                adjList = main_graph.subgraph(molecule).copy()
                baseoil_molecules.append(adjList)
        break           
    
    carbon_number = {}
    for i, base in enumerate(baseoil_molecules):
        # remove Hydrogen
        carbon_skeleton = base.copy()
        nodes_to_remove=[node for node in carbon_skeleton.nodes if atypes[node]==1]
        carbon_skeleton.remove_nodes_from(nodes_to_remove)
        
        if i==0:
            src,long = man.longest_chain(carbon_skeleton)
            G1 = carbon_skeleton.copy()
            carbon_number|=man.number_nodes_bfs(carbon_skeleton, src[0])
        else:
            corr_src = man.find_corresponding_node(G1, carbon_skeleton, src[0])
            carbon_number|=man.number_nodes_bfs(carbon_skeleton, corr_src)
    
    
    for baseoil_index in range(len(baseoil_molecules)):
        # baseoil_index = 2
        carbon_skeleton = baseoil_molecules[baseoil_index].copy()
        nodes_to_remove=[node for node in carbon_skeleton.nodes if atypes[node]==1]
        carbon_skeleton.remove_nodes_from(nodes_to_remove)
        
        # enu = sorted(['6-5', '8-7', '5-4', '4-3', '7-6', '12-11', '18-15', '19-16', '16-14',
        #        '11-10', '20-18', '17-14', '13-11', '15-13', '14-12', '30-29', '23-21',
        #        '25-23', '21-19', '9-8', '22-20', '10-9', '27-25', '29-27', '3-2',
               # '28-26', '26-24', '24-22', '2-1'])
        
        enu = ['2-1', '3-2', '4-3', '5-4', '6-5', '7-6', '8-7', '9-8', '10-9', '11-10', '12-11', '14-12', '16-14', '19-16', '21-19', '23-21', '25-23', '27-25', '13-11', '15-13', '18-15', '17-14', '20-18', '22-20', '24-22', '26-24', '28-26', '29-27', '30-29']
        
        data_buffer = []
        for step in bondorders:
            row_data = {}
            for atom1 in bondorders[step]:
                if atom1 not in carbon_skeleton.nodes:
                    continue
                for atom2 in bondorders[step][atom1]:
                    if atom2 not in carbon_skeleton.nodes:
                        continue
                    if atom1>atom2:
                        bond_order = bondorders[step][atom1][atom2]
                        a1 = carbon_number[atom1]
                        a2 = carbon_number[atom2]
                        column = f"{a1}-{a2}"
                        if column in enu:
                            num = enu.index(column)+1
                            row_data[num] = bond_order
            data_buffer.append(row_data)
            
        CC_bonds = pd.DataFrame(data_buffer).fillna(0)
        CC_bonds.index = bondorders.keys() 
        
        timestep = 0.25
        skipts   = (1600-300)/4
        n_xticks = 4
        
        species = CC_bonds.T.copy().iloc[:,:]   # copy is important
        species.columns = species.columns*timestep/1000-skipts
        species    = species.loc[:,0:] # skip the ramping steps
        
        ## skip rows stated with zero
        species = species[species.iloc[:,0] != 0]
        species = species.sort_index(axis=0)
        
        s = species.sum(axis=1)
        for idx, item in zip(s.index,s.values):
            if item<3500:
                counts[idx]+=1
        
        # xticks_freq = int(species.columns.size/n_xticks)
        # fig, ax = plt.subplots(figsize=(6,8))
        # sns.heatmap(species,cmap='jet',xticklabels=xticks_freq, ax=ax,
        #                                 vmin=0.0,vmax=2.0)
        
        # xticklabels = ax.get_xticklabels()
        # xticks = np.array([float(label.get_text()) for label in xticklabels])
        # ax.set_xticklabels(xticks.astype(int))
        # ax.set_yticklabels(ax.get_yticklabels(),rotation=0)
        
        # ax.set_xlabel('Time (ps)')
        # plt.savefig(f'heatmap_of_bonds_{baseoil_index}.png',dpi=500,bbox_inches='tight')
        # plt.show()

        # print(baseoil_index,counts)
#%%
plt.bar(counts.keys(),counts.values())
plt.yticks([0,5,10,15,20])
plt.xticks([0,5,10,15,20,25,30],color='red')
plt.xlabel("Bond Index",color='red')
plt.ylabel("Count")
plt.savefig('barplot.png',dpi=300)