# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 25 15:41:37 2024
"""

import MDAnalysis as mda
import numpy as np
from scipy.spatial.transform import Rotation as R

# Load the PDB file
input_pdb = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Material Structures\BHT.pdb"
output_xyz = "fake_trajectory.xyz"
u = mda.Universe(input_pdb)

# Calculate the center of mass
com = u.atoms.center_of_mass()

# Number of frames in the trajectory
num_frames = 360

# Open the trajectory file for writing
with open(output_xyz, 'w') as xyz_file:
    for frame in range(num_frames):
        # Translate molecule to the origin (for rotation)
        u.atoms.positions -= com
        
        # Define a rotation (e.g., incrementally rotate around the z-axis)
        angle = 360 / num_frames  # Rotate in equal increments
        rotation = R.from_euler('y', angle, degrees=True).as_matrix()
        
        # Apply the rotation
        u.atoms.positions = np.dot(u.atoms.positions, rotation.T)
        
        # Translate molecule back to the original position
        u.atoms.positions += com
        
        # Write the frame to the XYZ file
        xyz_file.write(f"{len(u.atoms)}\n")
        xyz_file.write(f"Frame {frame + 1}\n")
        for atom in u.atoms:
            xyz_file.write(f"{atom.type} {atom.position[0]:.3f} {atom.position[1]:.3f} {atom.position[2]:.3f}\n")

print(f"Fake trajectory saved to {output_xyz}")
