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
cutoff = 0.45
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff, mtypes=True,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

number_of_baseoils = 14
savedir = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway\pathway_BaseOils\{}".format(baseoil)
for base_oil_id in range(14,number_of_baseoils+1):
    img_counter = 1
    print(base_oil_id)
    outfilepath = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\tracker_{}_{}.txt'.format(baseoil,base_oil_id)
    smiles_list = []
    with open(outfilepath,'w') as outfile:
        print(cutoff)
        outfile.write(f"Bond Order Cutoff {cutoff}\n")
        last_frame = -1
        for frame,(step, neigh) in enumerate(neighbors.items()):
            main_graph = nx.Graph(neigh)
            molecules = bfp.get_molecules(neigh)
            carbon_count = 0
            for molecule in molecules:
                molecule_types = [mtypes[x] for x in molecule]
                carbon_count+=molecule_types.count(2)
                if base_oil_id in molecule_types:
                    
                    species = bfp.get_molecular_formula(molecule, atypes,
                                                        atomsymbols)
                    if species=="H62C30":
                        continue
                    
                    
                    moleculeGraph = main_graph.subgraph(molecule)
                    smiles = bfp.moleculeGraph2smiles(moleculeGraph, [1,6,8], None,
                                                      bo_analysis=False,atom_types=atypes)
                    
                    # if smiles in smiles_list or "." in smiles:
                    #     continue
                    
                    ## ploting
                    latex = bfp.make_molecular_formula_latex(species,sort=True)
                    title = "baseoil: {}, id: {}, frame: {},Count: {}\nSpecies={}".format(baseoil,base_oil_id,frame,img_counter,latex)
                    fig, ax = plt.subplots()
                    # man.draw_rdkit2D_from_graph(moleculeGraph,atypes,
                    #                             atomsymbols, ax=ax)
                    man.draw_rdkit2D_from_smiles(smiles,ax=ax)
                    ax.set_title(title)
                    outimg = "\\img_{}_id{}_{}".format(baseoil,base_oil_id,img_counter)
                    img_counter+=1
                    fig.savefig(savedir+outimg,dpi=300,bbox_inches='tight')
                    ############
                    
                    smiles_list.append(smiles)
                    if last_frame==frame:
                        line = str(frame)+": "+species+"\n"+smiles+"\n"
                        outfile.write(line)
                    else:
                        outfile.write("========================================\n")
                        line = str(frame)+": "+species+"\n"+smiles+"\n"
                        outfile.write(line)
                    last_frame = frame
            if carbon_count!=92:
                print('------------------------------------------')