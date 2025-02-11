# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Jan 28 03:37:49 2024
"""

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

baseoil = "PAO4"
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
cutoff = 0.35
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

n=25
skipts = (1600-300)/4
tick_fontsize  = 8
label_fontsize = 15
vmin   = 0
vmax   = 2
timestep = 0.25

fig, ax = plt.subplots(1,n+1,figsize=(2.25+n, 4),
                       gridspec_kw={'width_ratios': [1]*n+[0.1]})

flag = True
bonds_df = []

for i in range(n):
    firststep_neigh     = list(neighbors.values())[0]
    firststep_graph     = nx.Graph(firststep_neigh)
    molecule            = [k for k,v in mtypes.items() if v==i+1]
    
    subgraph = firststep_graph.subgraph(molecule)
    H_nodes  = [x for x in molecule if atypes[x]==1]
    backbone = subgraph.copy()
    backbone.remove_nodes_from(H_nodes)
    if flag:
        bonds = backbone.edges
        first_bonds = bonds
        flag = False
    else:
        bonds = []
        for a,b in first_bonds:
            bonds.append((a+92*i,b+92*i))
        current = backbone.edges
        print(i+1,set(bonds)==set(current))
        
    df = bfp.bondorder_evolution(bondorders, bonds, ts=timestep,
                                 skipts=skipts,plot='no')
    bonds_df.append(df)
    if i<n-1:
        hm = sns.heatmap(df,cmap='jet',xticklabels=2000,
                         ax=ax[i],vmin=0.0,vmax=2.0,cbar=False)
    else:
        cbar_ticks = np.arange(vmin,vmax+1,step=1,dtype=int)
        cbar_kws = {'ticks': cbar_ticks}
        sns.heatmap(df,cmap='jet',cbar_kws=cbar_kws,
                    xticklabels=2000,ax=ax[i],vmin=vmin,vmax=vmax,
                    cbar_ax=ax[-1])
        ax[-1].set_ylabel("Bond order", fontsize=label_fontsize)
        
    xticklabels = [int(float(x.get_text())) for x in ax[i].get_xticklabels()]
    ax[i].set_xticklabels(xticklabels,rotation=90,fontsize=tick_fontsize)
    yticks = [1]+list(range(5,df.index.size,5))+[df.index[-1]]
    ax[i].set_yticks(yticks)
    ax[i].tick_params(bottom=True, left=True, right=False, top=False)
    if i==0: ax[i].set_yticklabels(yticks,rotation=0,fontsize=tick_fontsize)
    # ax[i].set_title('{}A: {}'.format(i+1,subtitle[i]),
    #                 fontsize=label_fontsize-1)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution\bo_evolution_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
#%%---first bond breaking----
for df in bonds_df:
    for row in df.index:
        change = df.loc[row]<cutoff
        broke = change[change].index
        print(broke)
        break
    break