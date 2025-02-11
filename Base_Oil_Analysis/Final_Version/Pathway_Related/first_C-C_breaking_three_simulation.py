# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Jan 28 20:00:46 2024
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
simdir = ['Sim-1','Sim-2','Sim-3']

bonddata = {}
cutoff   = 0.35
for sim in simdir:
    directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(baseoil,baseoil,sim)
    bondfilepath = directory+'\\bonds.reaxc'
    bonddata[sim] = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%% mapping edges of all base oil to reference base oil (C-C bonds only)
ref_flag = True
bonds_df = {}

# loop through three simulations
for sim in simdir:
    neighbors = bonddata[sim]['neighbours']
    atypes    = bonddata[sim]['atypes']
    mtypes    = bonddata[sim]['mtypes']
    bondorders= bonddata[sim]['bondorders']
    atomsymbols = 'HCO'
    
    bonds_df[sim]=[]
    
    n=25
    skipts = (1600-300)/4
    tick_fontsize  = 8
    label_fontsize = 15
    vmin   = 0
    vmax   = 2
    timestep = 0.25


    # loop through 25 molecules
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
        center = nx.center(backbone)
        #getting reference graph
        if ref_flag:
            ref_id = i+1
            ref = backbone.copy()
            ref_flag = False
            
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
        bonds_df[sim].append(df)
        
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
            # color the center
            nx.draw_networkx_nodes(backbone, pos=pos, node_size=node_size,
                                   node_color='tab:grey', edgecolors='black',
                                   nodelist=center)
                    
            plt.text(0.4,-1.0,'$L_{8}$',fontsize=fontsize+15)
            plt.text(-0.9,-0.7,'$L_{11}$',fontsize=fontsize+15)
            plt.text(0.3,1.3,'$L_{10}$',fontsize=fontsize+15)
            rand = random.randint(0,1000000)
            savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs'
            plt.savefig(savedir+f'\\marking_edge_PAO4_{rand}.png', dpi=500,
                        bbox_inches='tight')
            plt.show()
#%%---first bond breaking----    
label_fontsize = 15
tick_fontsize = 15
print(f'Bond Order Cutoff: {cutoff}')
first_break = {x:0 for x in range(1,29+1)}

total_break = 0
branched_break = 0

for sim in simdir:
    print(sim)
    for i,df in enumerate(bonds_df[sim]):
        bond_break_times = []
        for mol_id in df.index:
            bond = df.loc[mol_id] # panda series
            
            # strict condition
            change = bond<cutoff
            broke = change[change].index
            if not broke.empty:
                bond_break_times.append(broke.min())
            else:
                bond_break_times.append(np.inf)
        minn = min(bond_break_times)
        bond_id = bond_break_times.index(minn)+1
        print(f'Mol-{i+1}: First break of bond {bond_id} at {minn}')        
        if minn!=np.inf:
            first_break[bond_id]+=1
            total_break+=1
            if 10<=bond_id<=15:
                branched_break+=1

plt.style.use('default')
fig, ax = plt.subplots(figsize=[8,6])
ax.bar(first_break.keys(),first_break.values(),color='#FFD54F',
        edgecolor='black')
ax.set_xlabel('Bond index',fontsize=label_fontsize)
ax.set_ylabel('C-C bond dissociation count',fontsize=label_fontsize)
ax.tick_params(axis='both', labelsize=tick_fontsize)

rand = random.randint(0,1000000)
savedir2 = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'
fig.savefig(savedir2+'\\bond_disso_count_{}.png'.format(rand),dpi=500,bbox_inches='tight')