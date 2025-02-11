# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Nov 29 04:09:53 2023
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

base_oil_id = 1
print("Base oil ID:",base_oil_id)
print("Bond order cutoff:",cutoff)

intermediates = {}
for frame,(step, neigh) in enumerate(neighbors.items()):
    intermediates[frame]=[]
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        molecule_types = [mtypes[x] for x in molecule]
        if base_oil_id in molecule_types:
            intermediates[frame].append(frozenset(molecule))
        
temp_life = {} 
for frame, molecules in intermediates.items():
    for molecule in molecules:
        if molecule in temp_life.keys():
            temp_life[molecule].append(frame)
        else:
            temp_life[molecule]=[frame]
#%%
life = {}
cutoff_life = 50
count=0
with open("life.txt","w") as lf:
    lf.write(f"# Base oil ID: {base_oil_id}, Bond order cutoff: {cutoff}\n")
    lf.write("species"+","+"life"+"\n")
    for molecule, lives in temp_life.items():
        species=bfp.get_molecular_formula(molecule, atypes, atomsymbols)
        lifetime = max(lives)-min(lives)
        if lifetime < cutoff_life:
            continue
        life[molecule]=lifetime
        lf.write(species+","+str(lifetime)+"\n")
        count+=1
        if species=='H2CO':
            print(species,"--",lifetime)
    lf.write(f"# Total count:{count}")
    print(f"Count {count}")