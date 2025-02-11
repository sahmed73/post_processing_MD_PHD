# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Sep  6 07:39:56 2024
"""

import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import numpy as np


topo = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0001\Production\Sim-1\onset.data'
trej = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0001\Production\Sim-1\onset.lammpstrj'
u = mda.Universe(topo, trej, format='LAMMPSDUMP')
#%%
for ts in u.trajectory:
    wrapped_coords = u.atoms.positions
    unwrapped_coords = np.copy(wrapped_coords)
    
    for i in range(1, len(wrapped_coords)):
        delta = wrapped_coords[i] - wrapped_coords[i-1]
        # Apply the PBC by removing jumps over box boundaries
        box_length = ts.dimensions[:3]
        delta = delta - np.rint(delta / box_length) * box_length
        unwrapped_coords[i] = unwrapped_coords[i-1] + delta
        
#%%
# File paths
topo = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0001\Production\Sim-1\onset.data'
trej = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0001\Production\Sim-1\onset.lammpstrj'

# Load universe
u = mda.Universe(topo, trej, format='LAMMPSDUMP')

# Output trajectory file name
output_trajectory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0001\Production\Sim-1\unwrapped_trajectory.dump'

# Open a writer to save the new trajectory
with mda.Writer(output_trajectory, n_atoms=u.atoms.n_atoms, format='LAMMPS') as W:
    for ts in u.trajectory:
        wrapped_coords = u.atoms.positions
        unwrapped_coords = np.copy(wrapped_coords)
        
        # Unwrap the coordinates manually
        for i in range(1, len(wrapped_coords)):
            delta = wrapped_coords[i] - wrapped_coords[i - 1]
            # Apply periodic boundary conditions (PBC)
            box_length = ts.dimensions[:3]  # Box dimensions
            delta = delta - np.rint(delta / box_length) * box_length
            unwrapped_coords[i] = unwrapped_coords[i - 1] + delta

        # Update positions in the universe with unwrapped coordinates
        u.atoms.positions = unwrapped_coords

        # Write the current frame to the new trajectory file
        W.write(u.atoms)

print(f"Unwrapped trajectory saved as: {output_trajectory}")

        
    
    





#%%
def unwrap_coordinates(wrapped_coords, box_length):
    """
    Unwraps the coordinates by removing periodic boundary effects.
    
    Parameters:
    - wrapped_coords: numpy array of shape (n_steps, n_atoms, 3) containing the wrapped coordinates.
    - box_length: The length of the simulation box (assuming cubic box).

    Returns:
    - unwrapped_coords: numpy array of the same shape as wrapped_coords with unwrapped coordinates.
    """
    unwrapped_coords = np.copy(wrapped_coords)
    
    # Iterate over time steps
    for i in range(1, len(wrapped_coords)):
        delta = wrapped_coords[i] - wrapped_coords[i-1]
        # Apply the periodic boundary conditions by removing jumps over box boundaries
        delta = delta - np.rint(delta / box_length) * box_length
        unwrapped_coords[i] = unwrapped_coords[i-1] + delta

    return unwrapped_coords

def load_wrapped_coordinates(file_path):
    """
    Loads the wrapped coordinates from a text file.
    Assumes file contains rows with x, y, z coordinates for each atom at each time step.

    Parameters:
    - file_path: Path to the file containing wrapped coordinates.

    Returns:
    - wrapped_coords: A numpy array of shape (n_steps, n_atoms, 3).
    """
    # Load coordinates from file
    return np.loadtxt(file_path).reshape(-1, n_atoms, 3)

def save_unwrapped_coordinates(unwrapped_coords, output_file):
    """
    Saves the unwrapped coordinates to a file.

    Parameters:
    - unwrapped_coords: A numpy array of unwrapped coordinates.
    - output_file: Path to the output file.
    """
    np.savetxt(output_file, unwrapped_coords.reshape(-1, 3), fmt='%.6f')


# Parameters
input_file = 'wrapped_coordinates.txt'  # Path to the input file with wrapped coordinates
output_file = 'unwrapped_coordinates.txt'  # Path to save the unwrapped coordinates
box_length = 10.0  # Length of the cubic box (modify according to your system)
n_atoms = 100  # Number of atoms (adjust based on your system)

# Load wrapped coordinates
wrapped_coords = load_wrapped_coordinates(input_file)

# Unwrap the coordinates
unwrapped_coords = unwrap_coordinates(wrapped_coords, box_length)

# Save the unwrapped coordinates
save_unwrapped_coordinates(unwrapped_coords, output_file)

print(f"Unwrapped coordinates saved to {output_file}")
