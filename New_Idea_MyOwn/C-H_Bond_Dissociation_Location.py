# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Mar  1 02:51:27 2025
"""

import magnolia.bondfile_parser as bfp


dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"
filename = r"\bonds.reaxc"

atom_connectivity = bfp.parsebondfile(dirr+filename, mtypes=True)
#%%
import matplotlib.pyplot as plt
import numpy as np

# Load connectivity and atomic type data
neighbors = atom_connectivity['neighbours']  # Adjacency list for each step
atypes    = atom_connectivity['atypes']      # Maps atom ID to atomic type (H=1, C=2, O=3, etc.)

# Step 1: Classify Carbons at Step 0
carbon_types = {}  # Store classification of each carbon
initial_graph = neighbors[0]  # Get the first step connectivity

for parent, children in initial_graph.items():
    if atypes[parent] == 2:  # Only classify Carbon (C=2)
        bonded_carbons = [x for x in children if atypes[x] == 2]  # Count bonded carbons
        num_carbons = len(bonded_carbons)

        if num_carbons == 1:
            carbon_types[parent] = "Primary"
        elif num_carbons == 2:
            carbon_types[parent] = "Secondary"
        elif num_carbons == 3:
            carbon_types[parent] = "Tertiary"

N_primary = list(carbon_types.values()).count('Primary')
N_secondary = list(carbon_types.values()).count('Secondary')
N_tertiary = list(carbon_types.values()).count('Tertiary')
count = dict()
count['Primary'] = N_primary
count['Secondary'] = N_secondary
count['Tertiary'] = N_tertiary

# Step 2: Count Hydrogens per Carbon at Any Step
def count_hydrogens(neigh, atypes):
    """Counts hydrogen atoms bonded to each carbon."""
    hydrogen_counts = {}

    for atom, bonded_atoms in neigh.items():
        if atypes[atom] == 2:  # Only count hydrogens around Carbon
            hydrogen_count = sum(1 for neighbor in bonded_atoms if atypes[neighbor] == 1)  # Count H
            hydrogen_counts[atom] = hydrogen_count

    return hydrogen_counts

# Initialize hydrogen count at step 0
previous_hydrogens = count_hydrogens(neighbors[0], atypes)  

# Step 3: Track Incremental Hydrogen Loss Over Time
H_loss_trend = {"Primary": [], "Secondary": [], "Tertiary": []}

for step, current_graph in neighbors.items():  # Loop through steps
    current_hydrogens = count_hydrogens(current_graph, atypes)
    
    # Track total H loss per carbon type at this step
    step_loss = {"Primary": 0, "Secondary": 0, "Tertiary": 0}

    for atom in current_hydrogens:
        if atom in previous_hydrogens:  # Ensure atom was present in the previous step
            h_lost = previous_hydrogens[atom] - current_hydrogens[atom]  # Compare with the last step

            if h_lost > 0:  # If hydrogens were lost
                c_type = carbon_types.get(atom, None)
                if c_type:
                    step_loss[c_type] += h_lost  # Sum H loss per step per type

    # Append total H loss per step
    for c_type in H_loss_trend:
        H_loss_trend[c_type].append(step_loss[c_type])
    
    # Update previous_hydrogens to current for the next iteration
    # previous_hydrogens = current_hydrogens  

# Step 4: Plot Incremental Hydrogen Loss Over Time 
steps = np.array(list(neighbors.keys()))
time = steps*0.25/1000

start = 325
end = None

fig, ax = plt.subplots(dpi=350)

for key, value in H_loss_trend.items():
    value = np.array(value)
    x     = time[time>325] - time[time>325].min()
    y     = value[time>325]
    ax.plot(x,y,label=key)

plt.xlabel("Time (ps)")
plt.ylabel("H Atoms Loss")
plt.legend()
plt.title("Incremental H Loss from Different\nCarbon Types Over Time\n", fontsize=12)
plt.show()

