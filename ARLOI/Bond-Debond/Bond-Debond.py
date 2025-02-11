# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 22 00:32:37 2024
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt

dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0001\Production\300-500K_TRR=1Kpps\Sim-1'
bondfile=dirr+r'\bonds.out'
atominfo=bfp.parsebondfile(bondfile,mtypes=True)
#%%
neighbours=atominfo['neighbours']
atypes=atominfo['atypes']
mtypes=atominfo['mtypes']
timestep=0.25

for track in range(1,20+1):
# track=1 # tacking molecule 1
    tracker={}
    print(f'--------------{track}-----------------')
    for step, neigh in neighbours.items():
        time = step*timestep/1000
        tracker[time]=[]
        molecules=bfp.get_molecules(neigh)
        for molecule in molecules:
            mmm = [mtypes[atom] for atom in molecule]
            most_freq=max(mmm, key=mmm.count)
            if most_freq==track:
                tracker[time].append(molecule)
    
    last = ""
    start = 0
    
    for time, molecules in tracker.items():
        if len(molecules)==0: continue
        molecule = list(molecules)[0]
        species = bfp.get_molecular_formula(molecule, atypes, 'HCO')
        
        if molecule == last:
            end = time
        else:
            if last:  # Avoid printing for the first iteration
                last_species=bfp.get_molecular_formula(last, atypes, 'HCO')
                print(f"{start}-{end}: {last_species}")
            start = time
            end = time  # Update end when species changes
        last = molecule
    
    # Print the last range after the loop
    if last:
        last_species=bfp.get_molecular_formula(last, atypes, 'HCO')
        print(f"{start}-{end}: {last_species}")

        





