# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 17 11:35:47 2023
"""

from collections import deque
import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt

directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"
bondfilepath = directory+'\\bonds.reaxc'
bonddata = bfp.parsebondfile(bondfilepath,mtypes=True)
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
#%%
from rdkit import Chem
from rdkit.Chem import Draw
def find_distances(adjacency_list, source_atom, atypes):
    """
    Finds the distance of each atom from the source atom using BFS algorithm.
    
    :param adjacency_list: A dictionary representing the graph (molecule), 
                           where keys are atoms and values are lists of adjacent atoms.
    :param source_atom: The atom from which distances are calculated.
    :return: A dictionary with each atom and its distance from the source atom.
    """
    distances = {atom: None for atom in adjacency_list}  # Initialize distances
    distances[source_atom] = 1  # Distance to source atom is 1

    queue = deque([source_atom])  # Queue for BFS

    while queue:
        current_atom = queue.popleft()

        for neighbor in adjacency_list[current_atom]:
            if distances[neighbor] is None:  # If not visited
                distances[neighbor] = distances[current_atom] + 1
                queue.append(neighbor)
                
    return distances


## selecting a proper source for each base oil molecule
firstneigh = neighbors[0]
proper_source_dict = {}
ref_flag = True
baseoil_molecule_list = [] # just for fun
for i in range(25):
    baseoil_molecule = {}
    for parent, children in firstneigh.items():
        # only considering carbon atoms as parent and children
        if mtypes[parent]==i+1 and atypes[parent]==2:
            carbon_children = [x for x in children if atypes[x]==2]
            baseoil_molecule[parent]=carbon_children
    
    
    max_dist = -1
    for source in baseoil_molecule:
        dist = find_distances(baseoil_molecule,source,atypes)
        if max(dist.values())>max_dist:
            max_dist = max(dist.values())
            if ref_flag:
                # ref list will be fixed for the first molecule
                ref = sorted(dist.values())
            
    for source in baseoil_molecule:
        dist = find_distances(baseoil_molecule,source,atypes)
        if max(dist.values())==max_dist and ref==sorted(dist.values()):
            proper_source = source
    baseoil_molecule_list.append(baseoil_molecule)
    
    proper_source_dict[i+1]=proper_source
    # print('Molecule {}: source={}, max dist={}'.format(i+1,proper_source,max_dist))
    ref_flag = False


savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\graphs'
figsize = [10,6]
node_size = 800
for molecule_ID in range(1,26):
    source = proper_source_dict[molecule_ID]
    
    firstStep_baseoil_molecule = {}
    for parent, children in firstneigh.items():
        # only considering carbon atoms as parent and children
        if mtypes[parent]==molecule_ID and atypes[parent]==2:
            carbon_children = [x for x in children if atypes[x]==2]
            firstStep_baseoil_molecule[parent]=carbon_children
    firstStep_baseoil_moleculeGraph = nx.Graph(firstStep_baseoil_molecule)
    dist = find_distances(firstStep_baseoil_molecule, source, atypes)
    
    for step,neigh in neighbors.items():
        baseoil_molecule = {}
        for parent, children in neigh.items():
            if mtypes[parent]==molecule_ID and atypes[parent]==2:
                carbon_children = [x for x in children if atypes[x]==2]
                baseoil_molecule[parent]=carbon_children
        baseoil_moleculeGraph = nx.Graph(baseoil_molecule)
        firstStep_edge_len = len(firstStep_baseoil_moleculeGraph.edges)
        edge_len = len(baseoil_moleculeGraph.edges)
        if not nx.is_isomorphic(firstStep_baseoil_moleculeGraph,baseoil_moleculeGraph) and firstStep_edge_len>edge_len:
            print("Molecule {}, Time {} ps".format(molecule_ID,step*0.25/1000))
            # color_map = ['r' if node==source else 'tab:blue' for node in firstStep_baseoil_moleculeGraph.nodes()]
            # plt.figure(figsize=figsize)
            # pos = nx.spring_layout(firstStep_baseoil_moleculeGraph)
            # nx.draw(firstStep_baseoil_moleculeGraph,pos=pos,
            #         node_color=color_map,with_labels=True, node_size=node_size)
            # plt.title("Molecule {}: Before".format(molecule_ID))
            # plt.savefig(savedir+'\\molecule_{}_Before.png'.format(molecule_ID),
            #             dpi=500, bbox_inches='tight')
            # plt.show()
            
            
            # plt.figure(figsize=figsize)
            # components = nx.connected_components(baseoil_moleculeGraph)
            # shift = [0.1,0.4]
            # shift_flag = True
            # for comp in components:
            #     print(len(comp))
            #     if shift_flag:
            #         pos1 = {x:y+shift for x,y in pos.items() if x  in comp}
            #         shift_flag=False
            #     else:
            #         pos1 = {x:y for x,y in pos.items() if x  in comp}
                    
            #     subgraph = baseoil_moleculeGraph.subgraph(comp)
            #     color_map = ['r' if node==source else 'tab:blue' for node in subgraph.nodes()]
            #     nx.draw(subgraph,pos=pos1,node_color=color_map,
            #             with_labels=True, node_size=node_size)
            
            
            # plt.title("Molecule {}: After".format(molecule_ID))
            # plt.savefig(savedir+'\\molecule_{}_After.png'.format(molecule_ID),
            #             dpi=500, bbox_inches='tight')
            # plt.show()
            
            smiles = bfp.moleculeGraph2smiles(firstStep_baseoil_moleculeGraph, [1,6,8],2, bo_analysis=False,
                                              atom_types=atypes.copy())
            molecule = Chem.MolFromSmiles(smiles)

            # Compute 2D coordinates
            Chem.rdDepictor.Compute2DCoords(molecule)
            
            # Draw the molecule and display it using matplotlib
            fig, ax = plt.subplots()
            img = Draw.MolToImage(molecule)
            plt.imshow(img)
            plt.axis('off')  # No axis for the image
            plt.title("Molecule {}: Before".format(molecule_ID))
            plt.show()
            
            
            
            smiles = bfp.moleculeGraph2smiles(baseoil_moleculeGraph, [1,6,8],2, bo_analysis=False,
                                              atom_types=atypes.copy())
            molecule = Chem.MolFromSmiles(smiles)

            # Compute 2D coordinates
            Chem.rdDepictor.Compute2DCoords(molecule)
            
            # Draw the molecule and display it using matplotlib
            fig, ax = plt.subplots()
            img = Draw.MolToImage(molecule)
            plt.imshow(img)
            plt.axis('off')  # No axis for the image
            plt.title("Molecule {}: After".format(molecule_ID))
            plt.show()
            break
            
            
            
            

