# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 01:37:20 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import os
import random
import math
import sys
from rdkit import Chem

def check_bond(ssbo,graph,bondtype,dbol=(1.1,2.1),tbol=(2.3,math.inf)):
    '''
    Parameters
    ----------
    ssbo : Dictionary of dictionary
        ssbo[atom_1][atom_2]=bond_order.
    component : any iterable
        nodelist
    bondtype : string
        'double', 'tripple' or 'both'
    dbol : TYPE, optional
        DESCRIPTION. The default is (1.3,2.3).
        dbol = double bond order limit
    tbol : TYPE, optional
        DESCRIPTION. The default is (2.3,math.inf).
        tbol = tripple bond order limit

    Returns
    -------
    None.

    '''
    edge_list = graph.edges
    weighted_edge_list = []
    if bondtype == 'double':
        for e in edge_list:
            a,b = e
            low, high = dbol
            if low<=ssbo[a][b]<high:
                weighted_edge_list.append((a,b))
                
    return weighted_edge_list


baseoil = "Squalane"
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1'.format(baseoil,baseoil)
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

neighbours,atomtypes,bondorder = bfp.get_neighbours(bondfile,bo=True)
#%%-------

# ran = random.randint(10000,99999)
# savedir = 'python_outputs\\pathway\\pathway_snaps'+'-'+str(ran)
# os.mkdir(savedir)

def graph_to_mol(G, atypes,asyms):
    mol = Chem.RWMol()

    # Add atoms
    node_to_idx = {}
    for node in G.nodes():
        atom_symbol = asyms[atypes[node]-1]
        atom = Chem.Atom(atom_symbol)
        atom_idx = mol.AddAtom(atom)
        node_to_idx[node] = atom_idx

    # Add bonds
    for edge in G.edges():
        start, end = edge[0], edge[1]
        mol.AddBond(node_to_idx[start], node_to_idx[end])

    return mol.GetMol()


molecules = []
atomsymbols = ['H','C','O']
seek_atoms= {2541, 1781, 1780, 1749} #H2CO
print(bfp.get_molecular_formula(seek_atoms, atomtypes, atomsymbols))
checklist= []
g = []
g_chem = []
double_bond = []
triple_bond = []
timestep = list(neighbours.keys())
print(len(neighbours))
for step,neigh in neighbours.items():
    G = nx.Graph(neigh)    
    connected = nx.connected_components(G)
    for component in connected:
        '''if bfp.get_formula(component, atomtypes)=='H2CO':
            print(component)'''
        if component & seek_atoms:
            common = component & seek_atoms
            if component not in checklist:
                s = G.subgraph(component)
                species = bfp.get_molecular_formula(component, atomtypes,atomsymbols)
                mol = graph_to_mol(s,atomtypes,atomsymbols)
                img =Chem.Draw.MolToImage(mol)
                img.save("{}_image.png".format(species))
                # g.append((s,step))
                # g_chem.append(bfp.get_molecular_formula(component, atomtypes,atomsymbols))
                # double_bond.append(check_bond(bondorder[step], s, 'double'))
                # checklist.append(component)
    
#%%----------

print('Done Timestep')
for i,sub_step in enumerate(g):
    sub,step = sub_step
    new_nodes = {}
    atomsymbols = 'HCO'
    atomocolors = ['lightgray','lightskyblue','lightcoral']
    color_map = []
    for node in sub.nodes():
        at = atomtypes[node]
        new_nodes[node] = atomsymbols[at-1]
        color_map.append(atomocolors[at-1])
    
    
    pos = nx.layout.spectral_layout(sub)
    #pos = nx.circular_layout(sub)
    pos = nx.spring_layout(sub, pos=pos, iterations=500)
    
    nx.draw_networkx(
        sub,
        pos,
        node_color=color_map,
        with_labels=True,
        labels=new_nodes,
    )
    
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=double_bond[i],
        width=8,
        alpha=0.5,
        edge_color="yellow",
    )
    species = bfp.make_molecular_formula_latex([g_chem[i]],sort=True)
    frame = 'Frame_{}_{}'.format(i,g_chem[i])
    title = 'Frame: {}'.format(i)+'        '+species+'       '+'Timestep: '+str(step)
    plt.title(title)
    plt.savefig(savedir+'\\'+frame, dpi=300,bbox_inches='tight')
    plt.show()
