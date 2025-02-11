# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov 21 23:20:58 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import sys
import matplotlib.pyplot as plt
import magnolia.molecular_analysis as man

baseoil = "PAO4"
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
bonddata = bfp.parsebondfile(bondfilepath,mtypes=True,bo=True)
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
#%%
number_of_baseoils = 25
# base_oil_id = 2 # 1 to 25
savedir = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway\pathway_BaseOils"
for base_oil_id in range(1,2):
    img_counter = 0
    print(base_oil_id)
    outfilepath = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\tracker_{}_{}.txt'.format(baseoil,base_oil_id)
    smiles_list = []
    with open(outfilepath,'w') as outfile:
        last_frame = -1
        for frame,(step, neigh) in enumerate(neighbors.items()):
            main_graph = nx.Graph(neigh)
            molecules = bfp.get_molecules(neigh)
            for molecule in molecules:
                molecule_types = {mtypes[x] for x in molecule}
                if base_oil_id in molecule_types:
                    species = bfp.get_molecular_formula(molecule, atypes, 'HCO')
                    if species=="H62C30":
                        continue
                    moleculeGraph = main_graph.subgraph(molecule)
                    smiles = bfp.moleculeGraph2smiles(moleculeGraph, [1,6,8], None,
                                                      bo_analysis=False,atom_types=atypes)
                    if smiles in smiles_list or "." in smiles:
                        continue
                    # print(smiles)
                    title = "baseoil: {}, id: {}, frame: {}\nSpecies={}".format(baseoil,base_oil_id,frame, species)
                    fig, ax = plt.subplots()
                    man.draw_rdkit2D(smiles, title=title, ax=ax)
                    ax.set_title(title)
                    outimg = "\\img_{}_id{}_{}".format(baseoil,base_oil_id,img_counter)
                    img_counter+=1
                    fig.savefig(savedir+outimg,dpi=300,bbox_inches='tight')
                    smiles_list.append(smiles)
                    if last_frame==frame:
                        line = str(frame)+": "+species+"\n"+smiles+"\n"
                        outfile.write(line)
                    else:
                        outfile.write("========================================\n")
                        line = str(frame)+": "+species+"\n"+smiles+"\n"
                        outfile.write(line)
                    last_frame = frame