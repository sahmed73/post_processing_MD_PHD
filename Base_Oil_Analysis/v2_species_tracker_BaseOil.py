# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov  7 23:36:53 2023
"""

import magnolia.bondfile_parser as bfp
import sys
import matplotlib.pyplot as plt
import random
import numpy as np
import networkx as nx

### Directory ###
base_oil = 'PAO4'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(base_oil,base_oil)
filename  = "\\bonds.reaxc"
bondfilepath = directory+filename

### Parsing Bondfile ###
bonddata = bfp.parsebondfile(bondfilepath)
#%% Assigning Molecule Identification Number
neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
atomsymbols= ['H','C','O']

mol_to_molUID = {} # key=tuple(sorted(molecule)), value=molUID
molUID_to_mol = {} # key=molUID, value=molecule

molUID = 1 # molecule uniqueID
nodeID = 1 # node_ID
reaction_network = nx.DiGraph()
previous_molecules = []
for frame, (step, neigh) in enumerate(neighbours.items()):
    natoms             = len(neigh)
    atom_checklist     = {}
    molecules          = bfp.get_molecules(neigh)
    for molecule in molecules:
        if molecule in previous_molecules:
            key = tuple(sorted(molecule))
            pre_molUID = mol_to_molUID[key]
            reaction_network.add_node(nodeID,molUID=pre_molUID,frame=frame)
            nodeID+=1
        else:
            for pre_mol in previous_molecules:
                shared_atoms = pre_mol & molecule
                if shared_atoms:
                    key = tuple(sorted(pre_mol))
                    pre_molUID = mol_to_molUID[key]
                    
                    pass
                else:
                    pass
            
            