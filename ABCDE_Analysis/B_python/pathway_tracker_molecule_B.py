# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec 29 16:13:59 2023
"""

import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import time
import magnolia.needless_essential as ne
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem


### Directory ###
molecule = 'molecule_B'
# directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1000K\Sim-1"
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\Pure_50B\Production\1050\Sim-1'
filename  = "\\bonds.reaxc"
bondfile = directory+filename
timestep  = 0.25
cutoff = 0.3

bonddata = bfp.parsebondfile(bondfile, cutoff=cutoff, bo=True, mtypes=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
mtypes    = bonddata['mtypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

seek = "H2O"

seek_molecules = set()
for step, neigh in neighbors.items():
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        species = bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        if species==seek and frozenset(molecule) not in seek_molecules:
            seek_molecules.add(frozenset(molecule))
            print(molecule)
            for atom1 in molecule:
                for atom2 in molecule:
                    if atom1>atom2:
                        sym1 = atomsymbols[atypes[atom1]-1]
                        sym2 = atomsymbols[atypes[atom2]-1]
                        bb = bondorders[step][atom1].get(atom2,None)
                        if bb:
                            print((sym1,sym2),bb)
#%%
seek_molecules_copy = seek_molecules.copy()
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway'
with open(savedir+f"\\life_{seek}_{molecule}.txt","w") as lf:
    lf.write(f"# Pathway tracker for {seek} in "+
             f"{molecule} oxidation, Bond order cutoff: {cutoff}\n")
    lf.write("species"+","+"life"+"(starting_frame,ending_frame)"+
             ",atom_contributuin"+"\n")
    
    while seek_molecules_copy:
        seek_molecule = seek_molecules_copy.pop()
        print(seek_molecule)
        
        pathways = []
        for frame,(step, neigh) in enumerate(neighbors.items()):
            molecules = bfp.get_molecules(neigh)
            for molecule in molecules:
                common = molecule & seek_molecule
                chem = bfp.get_molecular_formula(common, atypes, atomsymbols)
                if common:
                    pathways.append((frame,frozenset(molecule)))
        
        temp_life = {} 
        for frame, molecule in pathways:
            if molecule in temp_life.keys():
                temp_life[molecule].append(frame)
            else:
                temp_life[molecule]=[frame]
        
        life = {}
        cutoff_life = 50
        count=1
        
        for molecule, lives in temp_life.items():
            molecule_ids = {}
            species=bfp.get_molecular_formula(molecule, atypes, atomsymbols)
            lifetime = max(lives)-min(lives)
            minn = min(lives)
            maxx = max(lives)
            if lifetime < cutoff_life:
                continue
            for atom in molecule:
                molecule_type = mtypes[atom]
                if molecule_type in molecule_ids:
                    molecule_ids[molecule_type].append(atom)
                else:
                    molecule_ids[molecule_type]=[atom]
                    
            for key, value in molecule_ids.items():
                value_species = bfp.get_molecular_formula(value, atypes, atomsymbols)
                molecule_ids[key]=value_species
            life[molecule]=lifetime
            lf.write(f"{count}. "+species+","+f" ({minn},{maxx},{lifetime}), {molecule_ids}"+"\n")
            count+=1
        lf.write(f"# Total count:{count}\n")
        lf.write("######\n\n")
        print(f"Count {count}")