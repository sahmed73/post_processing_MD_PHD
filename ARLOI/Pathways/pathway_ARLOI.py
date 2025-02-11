# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 18:33:10 2023

@author: Shihab
"""
import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import time
import magnolia.needless_essential as ne
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import magnolia.access_ucm_cluster as ucm



### Directory ###
base_oil = 'A000'
directory = "/mnt/borgstore/amartini/sahmed73/ARLOI-V1/4_ARLOI_24-09-07--22-34-03/20_PAO-OH_15_A0001/Production/Sim-1"
filename  = "/bonds.reaxc"
bondfile = ucm.local_copy(directory+filename)
timestep  = 0.25

atominfo = bfp.parsebondfile(bondfile, bo=True)

neighbours = atominfo['neighbours']
atomtypes  = atominfo['atypes']
bondorders = atominfo['bondorders']
#%%----------------------------
start_time = time.time()
savedir=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\expression_selectors'
exp_file = '{}\\exp-{}.txt'.format(savedir,base_oil)
bfp.expression_selector(neighbours, atomtypes,atomsybols='HCO',file=exp_file,heading=directory)
seek_molecule = [288,298,299,285]

ne.print_runtime(time.time()-start_time)
#%% run parsebondfile_asGraph
atomConnectivity = bfp.parsebondfile_asGraph(bondfile)
#%%
## draw_molecule using graph
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'
atomsymbols = 'HCO'
checklist = []
count=0
for frame, (step,neigh) in enumerate(neighbours.items()):
    molecules=bfp.get_molecules(neigh)
    for molecule in molecules:
        species = bfp.get_molecular_formula(molecule, atomtypes, atomsymbols)
        if species=='H60C30' and molecule not in checklist:
            atomConnectGraph = atomConnectivity[step]
            subgraph=atomConnectGraph.subgraph(molecule)
            plt.figure(figsize=[15,12])
            bfp.draw_molecule(subgraph,atomsymbols)
            plt.title('Frame {}'.format(frame),fontsize=20)
            plt.savefig(savedir+'\\draw_molecule_{}'.format(frame),
                        dpi=300, bbox_inches='tight')
            plt.show()
            for atom1 in molecule:
                for atom2 in molecule:
                    if bondorders[step][atom1].get(atom2,False) and atom1>atom2:
                        a = atomsymbols[atomtypes[atom1]-1]
                        b = atomsymbols[atomtypes[atom2]-1]
                        print("({},{})={}".format(a,b,
                                            bondorders[step][atom1][atom2]))
            print(frame,molecule)
            print("-------------------")
            checklist.append(molecule)
            count+=1
print(count)
#%% Draw molecule using RDKit
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway'
atomic_num = [1,6,8] # H,C,O
atomsymbols = 'HCO'
seek = 'H2O'
checklist = []
for frame, (step,neigh) in enumerate(neighbours.items()):
    molecules=bfp.get_molecules(neigh)
    for molecule in molecules:
        species = bfp.get_molecular_formula(molecule, atomtypes, atomsymbols)
        if species==seek and molecule not in checklist:
            atomConnectGraph = atomConnectivity[step]
            moleculeGraph=atomConnectGraph.subgraph(molecule)
            smiles = bfp.moleculeGraph2smiles(moleculeGraph, atomic_num,
                                              n_clusters=1, plot_cluster=True)
            print(smiles)
            print(molecule)
            checklist.append(molecule)
#%% Track a moleucle
checklist = []
atomic_num = [1,6,8] # H,C,O
atomsymbols = 'HCO'
track = {515, 2412, 2286, 2234, 2430} # H2CO2

track_species = bfp.get_molecular_formula(track, atomtypes, atomsymbols)

savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway\pathway_{}_{}_No_H'.format(base_oil,track_species)
if not os.path.exists(savedir):
    os.makedirs(savedir)

last_species = ""
count = 1
for frame, (step,neigh) in enumerate(neighbours.items()):
    molecules=bfp.get_molecules(neigh)
    for molecule in molecules:
        if molecule==track:
            continue
        shared = molecule & track
        species=bfp.get_molecular_formula(molecule, atomtypes, atomsymbols)
        # print(list(shared).count(2),list(shared).count(3))
        shared_types = [atomtypes[x] for x in shared]
        if shared_types.count(2)>0 or shared_types.count(3)>0:
            if molecule not in checklist and species!=last_species:
                atomConnectGraph = atomConnectivity[step]
                moleculeGraph = atomConnectGraph.subgraph(molecule)
                shared_part=bfp.get_molecular_formula(shared, atomtypes, atomsymbols)
                
                try:
                    smiles = bfp.moleculeGraph2smiles(moleculeGraph, atomic_num,
                                                      n_clusters=2)
                    mol = Chem.MolFromSmiles(smiles)
                    AllChem.Compute2DCoords(mol)
                except Exception as e:
                    print(f"An error occurred: {e}")
                    # If an error occurs, run the alternate function
                    smiles = bfp.moleculeGraph2smiles(moleculeGraph, atomic_num,
                                                      n_clusters=1)
                    mol = Chem.MolFromSmiles(smiles)
                    AllChem.Compute2DCoords(mol)
                    
                
                rdkit_img = Draw.MolToImage(mol,size=(1200, 1200))
               
                # Convert to Matplotlib figure
                fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
                ax.imshow(rdkit_img, aspect='equal')
                ax.axis('off')
                
                latex_text = bfp.make_molecular_formula_latex([species],sort=True)
                ax.text(0.5, 0.1, latex_text, fontsize=20, va='center', ha='center', transform=ax.transAxes)
                
                # Save the figure
                plt.savefig(savedir+'\\{}_{}_without_H_comparison.png'.format(count,species),
                            bbox_inches='tight', dpi=500)
                    
                count+=1

                checklist.append(molecule)
                last_species=species
            