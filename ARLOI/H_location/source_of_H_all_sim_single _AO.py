# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 11 03:12:12 2024

1. phenolic_OH
2. non_phenolic_OH
3. aromatic_H
4. aliphatic_H
5. self_exchange

"""

import magnolia.bondfile_parser as bfp
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

def source_of_H_count(atominfo):
    neighbours = atominfo['neighbours']
    atypes = atominfo['atypes']
    mtypes = atominfo['mtypes']

    all_H_sources = {key: [] for key in ['Ph—OH',
                                         'non-Ph—OH',
                                         'Ar',
                                         'Self',
                                         't-Bu',
                                         'Others']}
    
    scavenged_H_sources_count = dict.fromkeys(all_H_sources, 0)
    
    # Using the first step only for grouping
    first_neigh=list(neighbours.values())[0]
    G = nx.Graph(first_neigh)
    
    ## gorup the aromatic carbons
    cycles = [cycle for cycle in nx.cycle_basis(G) if len(cycle) == 6]
    C_aromatic = {node for cycle in cycles for node in cycle}
    
    ## group the tert-Butyl carbons
    C_tBu = set()
    for parent, children in first_neigh.items():
        if atypes[parent] == 2:
            # if all four children carbon
            if mtypes[parent]>20 and all(atypes[x]==2 for x in children) and len(children)==4:
                C_tBu.add(parent)
                for child in children:
                    if child not in C_aromatic: C_tBu.add(child)
    

    for parent, children in first_neigh.items():
        if atypes[parent] == 1:  # Check only H atoms
            child = children[0]  # assuming H only has one bond

            if mtypes[parent] <= 20: # if the molecule is PAO radical
                all_H_sources['Self'].append(parent)
            elif atypes[child] == 3:  # connected to oxygen
                second_children = first_neigh[child]
                if any(s_child in C_aromatic for s_child in second_children):
                    all_H_sources['Ph—OH'].append(parent)
                else:
                    all_H_sources['non-Ph—OH'].append(parent)
            elif child in C_aromatic:
                all_H_sources['Ar'].append(parent)
            elif child in C_tBu:
                all_H_sources['t-Bu'].append(parent)
            else:
                all_H_sources['Others'].append(parent)
                
    # looping for counting
    checked = set()
    for frame, (step, neigh) in enumerate(neighbours.items()):
        molecules = bfp.get_molecules(neigh)
        for molecule in molecules:
            skeleton_length = sum(1 for c in molecule if atypes[c] == 2)
            if skeleton_length == 30:  # Check only specified molecules
                for atom in molecule:
                    if mtypes[atom] <= 20 and atypes[atom] == 3:
                        for child in neigh[atom]:
                            if atypes[child] == 1 and child not in checked:
                                category = next(
                                    (key for key in all_H_sources if child in all_H_sources[key]), 'Others'
                                )
                                scavenged_H_sources_count[category] += 1
                                checked.add(child)

    return scavenged_H_sources_count

# File paths
all_sim_scavenged_H_sources_count={}
for sim in range(1,1+1):
    dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0003\Production\300-500K_TRR=1Kpps"
    filename = r"\bonds.out"
    bondfile_path = dirr + f'\\Sim-{sim}' + filename
    atominfo = bfp.parsebondfile(bondfile_path, mtypes=True)
    
    all_sim_scavenged_H_sources_count[f'Sim-{sim}'] = source_of_H_count(atominfo)
#%%    
df = pd.DataFrame(all_sim_scavenged_H_sources_count)
df.to_csv(f"{dirr}/scavenged_H_counts.csv")
# df=pd.read_csv(f"{dirr}/scavenged_H_counts.csv", index_col=0)
series=df.sum(axis='columns')

fig, ax = plt.subplots(dpi=350)
sns.barplot(x=series.index, y=series.values, ax=ax)

# Customize the plot using `ax`
ax.set_xticklabels(df.index, rotation=45)
ax.set_xlabel("H Source Category")
ax.set_ylabel("Count")
ax.set_title("Counts of Scavenged H by Source Category")
# fig.savefig(output_path, bbox_inches='tight', dpi=350)
