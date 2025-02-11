# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 30 20:50:35 2023
"""

## update species pathway tracker

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

seek = "H2CO2"

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
with open(savedir+f"\\life_{seek}.txt","w") as lf:
    lf.write(f"# Pathway tracker for {seek} in "+
             f"{baseoil} oxidation, Bond order cutoff: {cutoff}\n")
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