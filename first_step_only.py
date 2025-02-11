# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:18:00 2023

@author: arup2
"""

import numpy as np
import magnolia.bondfile_parser as bfp
import magnolia.dumpfile_parser as dfp
import sys
import random
import networkx as nx
import matplotlib.pyplot as plt

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K\Sim-1'
filename    = 'bonds.reaxc'
bondfile_path    = directory+'\\'+filename
dumpfile_path    = directory+'\\'+'oxidation.lammpstrj'

bonddata = bfp.parsebondfile(bondfile_path, bo=True, mtypes=True, fso=True)
dumpdata = dfp.parsedumpfile(dumpfile_path,fso=True)

neighbours  = bonddata['neighbours']
atomtypes   = bonddata['atypes']
atomsymbols = 'HCO'
bondorders  = bonddata['bondorders']
atompos     = dumpdata['position']

lignin = {'A':'H12C10O3',
           'B':'H34C26O4',
           'C':'H44C29O2',
           'D':'H30C19O3',
           'E':'H14C11O4'}

L = 'D'

plt.rcParams['font.family'] = 'Times New Roman'
count = {'O2':0,lignin[L]:0}
once = True
for step, neigh in neighbours.items():
    molecules = bfp.get_molecules(neigh)
    graph = nx.Graph(neigh)
    for molecule in molecules:
        formula = bfp.get_molecular_formula(molecule, atomtypes, atomsymbols)
        count[formula]+=1
        
        if formula==lignin[L] and once:
            # draw molecule
            subgraph = graph.subgraph(molecule)
            new_nodes = {}
            atomsymbols = 'HCO'
            atomocolors = ['lightgrey','black','red']#['#ffffff','#909090','#ff0d0d']
            atomsize    = [700,1800,1500]
            color_map = []
            node_size = []
            for node in subgraph.nodes():
                at = atomtypes[node]
                new_nodes[node] = atomsymbols[at-1]
                color_map.append(atomocolors[at-1])
                node_size.append(atomsize[at-1])
            
            benzene = nx.cycle_basis(subgraph)[0]
            benzene_graph = graph.subgraph(benzene)
            benzene_pos = nx.spring_layout(benzene_graph)
            
            pos = nx.circular_layout(subgraph)
            # pos.update(benzene_pos)
            pos = nx.spring_layout(subgraph, iterations=1200, pos=pos)
            
            plt.figure(figsize=(35,25))
            
            # pos.update(benzene_pos)
            
            draw_nodes = nx.draw_networkx_nodes(
                subgraph,
                pos,
                node_color=color_map,
                node_size=node_size,
            )
            draw_nodes.set_edgecolor('k')
                        
            nx.draw_networkx_edges(subgraph, pos,width=4)
            
            for edge in subgraph.edges:
                atom1, atom2 = edge
                bo = bondorders[0][atom1].get(atom2,None)
                
                if bo>1.2: 
                    nx.draw_networkx_edges(subgraph,pos,edgelist=[edge],
                                       edge_color='k',width=13)
                    nx.draw_networkx_edges(subgraph,pos,edgelist=[edge],
                                       edge_color='white',width=5)
                        
            latex = bfp.make_molecular_formula_latex([formula])
            plt.title(*latex,fontsize=55)
            plt.scatter([],[],label='H',color='#ffffff',s=700,edgecolors='k',linewidths=0.5)
            plt.scatter([],[],label='O',color='#ff0d0d',s=1500,edgecolors='k',linewidths=0.5)
            plt.scatter([],[],label='C',color='#909090',s=1800,edgecolors='k',linewidths=0.5)
            
            plt.legend(fontsize=50)
            
            plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs\graph_{}'.format(formula+'-'+str(random.randint(0,100000))), bbox_inches='tight',dpi=400)
            
            
            once = False

print(count)

