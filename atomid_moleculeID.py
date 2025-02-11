# -*- coding: utf-8 -*-
"""
Created on Wed May 10 00:04:42 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import magnolia.needless_essential as ne


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_O2_Mixture\AO1\50_AO1_200_O2 Model'

filename = r'\50_AO1_200_O2.data'

datafile = directory+filename


atom_symbols = 'HCO'

molecule_id_atom_id  = {}
atom_id_molecule_id  = {}
atom_id_type = {}
atom_id_symbols = {}

with open(datafile,'r') as f:
    for line in f:       
        atom_info = line.split()
        if len(atom_info)!=7:
            continue
        
        atom_id   = int(atom_info[0])
        mol_id    = int(atom_info[1])
        atom_type = int(atom_info[2])
        
        if mol_id not in molecule_id_atom_id.keys():
            molecule_id_atom_id[mol_id]=[]
            
        molecule_id_atom_id[mol_id].append(atom_id)
        atom_id_molecule_id[atom_id]=mol_id        
        atom_id_type[atom_id]=atom_type
        atom_id_symbols[atom_id]=atom_symbols[atom_type-1]

###########
seek_atom_in_molecule = {}

for atom_id,symbol in atom_id_symbols.items():
    mol_id = atom_id_molecule_id[atom_id]
    if symbol=='O' and len(molecule_id_atom_id[mol_id])==25:
        
        if mol_id not in seek_atom_in_molecule.keys():
            seek_atom_in_molecule[mol_id]=[]
        seek_atom_in_molecule[mol_id].append(atom_id)

oxygen_1 = []
oxygen_2 = []
oxygen_3 = []
for mol_id,atoms in seek_atom_in_molecule.items():
    oxygen_1.append(atoms[0])
    oxygen_2.append(atoms[1])
    oxygen_3.append(atoms[2])
    
print(ne.atoms2expression_selection(oxygen_1))













        