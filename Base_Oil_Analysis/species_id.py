# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Dec 20 14:55:52 2023
"""

import magnolia.bondfile_parser as bfp

baseoil = "PAO4"
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1".format(baseoil,baseoil)
bondfilepath = directory+'\\bonds.reaxc'
cutoff = 0.45
bonddata = bfp.parsebondfile(bondfilepath,cutoff=cutoff ,bo=True)
#%%
neighbors = bonddata['neighbours']
atypes    = bonddata['atypes']
bondorders= bonddata['bondorders']
atomsymbols = 'HCO'

timestep = 0.25
initial_temp = 300
final_temp = 1600
temp_ramprate = 4
skip_frame = int(((final_temp-initial_temp)/temp_ramprate)/timestep)

speciesID = {}
life_frame = {}

## asign a species id to each and every species
uniqueID = 1
once = True
for frame, (step, neigh) in enumerate(neighbors.items()):
    if frame <= skip_frame: continue
    actual_frame = frame - skip_frame
    
    molecules = bfp.get_molecules(neigh)
    for molecule in molecules:
        frozen_molecule = frozenset(molecule)
        if frozen_molecule not in speciesID:
            speciesID[frozen_molecule] = uniqueID
            uniqueID+=1
            life_frame[frozen_molecule] = (actual_frame,0)
        else:
            life_frame[frozen_molecule] = (life_frame[frozen_molecule][0],actual_frame)