# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Jan 26 01:50:25 2024
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man
import numpy as np
import seaborn as sns
import random
from collections import Counter


baseoil = "PAO4"
simdirr = ['Sim-1','Sim-2','Sim-3']
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
cutoff = 0.35
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%% mapping edges of all base oil to reference base oil (C-C bonds only) 
sim = 'Sim-1'
neighbors = bonddata[sim]['neighbours']
atypes    = bonddata[sim]['atypes']
mtypes    = bonddata[sim]['mtypes']
bondorders= bonddata[sim]['bondorders']
atomsymbols = 'HCO'


savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs'
n=25
skipts = (1600-300)/4
tick_fontsize  = 8
label_fontsize = 15
vmin   = 0
vmax   = 2
timestep = 0.25

flag = True
bonds_df = []

for i in range(n):
    firststep_neigh     = list(neighbors.values())[0]
    firststep_graph     = nx.Graph(firststep_neigh)
    # (i+1)th base oil molecule (1-25)
    molecule            = [k for k,v in mtypes.items() if v==i+1]
    
    # molecule graph
    subgraph = firststep_graph.subgraph(molecule)
    H_nodes  = [x for x in molecule if atypes[x]==1]
    backbone = subgraph.copy()
    # backbone of molecule-{i+1}
    backbone.remove_nodes_from(H_nodes)
    
    #getting reference graph
    if flag:
        ref_id = i+1
        ref = backbone.copy()
        flag = False
        
    # graph matcher object
    GM = nx.isomorphism.GraphMatcher(ref, backbone)
    if GM.is_isomorphic():
        mapped_edges = []
        mapping = GM.mapping
        for u,v in ref.edges():
            u_map,v_map = mapping[u], mapping[v]
            mapped_edge = (u_map,v_map)
            mapped_edges.append(mapped_edge)
            print(f'({u},{v})={mapped_edge}',end='\t')
        print(f'\n========={i+1}=========')
        
    else:
        raise ValueError(f'{i+1} is not isomorphic to {ref}')
        
        
    df = bfp.bondorder_evolution(bondorders, mapped_edges, ts=timestep,
                                 skipts=skipts,plot='no')
    bonds_df.append(df)
    
    ## Draw 1st graphs
    if i==0:
        node_color = []
        for node in backbone.nodes:
            if atypes[node]==3: node_color.append('red')
            else: node_color.append('tab:blue')
            
        pos = nx.kamada_kawai_layout(backbone)
        for key,value in pos.items():
            pos[key]=2*value
            
        plt.figure(figsize=[8,6])
        node_size = 150
        nx.draw(backbone,pos=pos,node_size=node_size,node_color=node_color)
        
        title = f'{i+1}'
        plt.title(title,fontsize=30)
        
        for k, (u, v) in enumerate(ref.edges):
            u_map,v_map = mapping[u], mapping[v]
            
            A = pos[u_map]
            B = pos[v_map]
            C = (A+B)/2 #+ np.array([-0.16,-0.06])
            
            L8 = [mapping[x] for x in range(23,30+1)]+[mapping[11]]
            L10 = [mapping[x] for x in range(1,10+1)]+[mapping[11]]
            L11 = [mapping[x] for x in range(12,22+1)]+[mapping[11]]
            
            
            fontsize = 12
            if u_map in L8 and v_map in L8:
                # C+=np.array([0.1,0.1])
                shift = [0.0,0.0]
                pos_L8 = C+np.array(shift)
                plt.text(*pos_L8,f'{k+1}',color='red',fontsize=fontsize)
                nx.draw_networkx_nodes(backbone, pos=pos, node_size=node_size,
                                       node_color='tab:red', edgecolors='black',
                                       nodelist=[u,v])
            elif u_map in L10 and v_map in L10:
                shift = [0.01,-0.06]
                pos_L10=C+np.array(shift)
                plt.text(*pos_L10,f'{k+1}',color='darkgreen',fontsize=fontsize)
                nx.draw_networkx_nodes(backbone, pos=pos, node_size=node_size,
                                       node_color='tab:green', edgecolors='black',
                                       nodelist=[u,v])
            elif u_map in L11 and v_map in L11:
                shift = [-0.05,0.04]
                if {u,v}=={13,22}: shift = [0,0]
                pos_L11 = C+np.array(shift)
                plt.text(*pos_L11,f'{k+1}',color='blue',fontsize=fontsize)
                nx.draw_networkx_nodes(backbone, pos=pos, node_size=node_size,
                                       node_color='tab:blue', edgecolors='black',
                                       nodelist=[u,v])
            else:
                print('something is wrong',(u_map,v_map))
                
        plt.text(0.4,-1.0,'$L_{8}$',fontsize=fontsize+15)
        plt.text(-0.9,-0.7,'$L_{11}$',fontsize=fontsize+15)
        plt.text(0.3,1.3,'$L_{10}$',fontsize=fontsize+15)
        plt.savefig(savedir+'\\marking_edge_PAO4.png', dpi=500,
                    bbox_inches='tight')
        plt.show()
#%%---first bond breaking----
print(f'Bond Order Cutoff: {cutoff}')
first_break = []
for i,df in enumerate(bonds_df):
    bond_break_times = []
    for mol_id in df.index:
        edge_list = df.loc[mol_id] # panda series
        
        # strict condition
        change = edge_list<cutoff
        broke = change[change].index
        if not broke.empty:
            bond_break_times.append(broke.min())
        else:
            bond_break_times.append(np.inf)
    minn = min(bond_break_times)
    num = bond_break_times.index(minn)+1
    print(f'Mol-{i+1}: First break of bond {num} at {minn}')        
    first_break.append(num)
    
bond_break_count = Counter(first_break)
plt.bar(bond_break_count.keys(),bond_break_count.values())
