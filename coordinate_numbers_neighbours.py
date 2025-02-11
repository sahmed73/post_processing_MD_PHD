# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Jan 26 14:07:07 2025
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Load the data file (e.g., XYZ or LAMMPS data)
u = mda.Universe(r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0003\DataFile\Single_Molecules\A0003_pre.xyz")

# Define the atom of interest (e.g., by name, type, or index)
atom_of_interest = u.select_atoms("name O")[0]  # Select the first oxygen atom
cutoff_radius = 3.0  # Define the cutoff radius in Å

# Select all neighboring atoms within the cutoff radius
neighbors = u.select_atoms(f"around {cutoff_radius} index {atom_of_interest.index}")

# Define van der Waals radii (in Å)
vdw_radii = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    # Add more elements as needed
}

# Define colors for each atom type
atom_colors = {
    "H": "gray",
    "C": "blue",
    "N": "green",
    "O": "red",
    "S": "yellow",
}

# Assign radii and colors to neighbors based on their element type
neighbor_radii = [
    vdw_radii.get(atom.element, 1.50)  # Use a default radius if element is not in the dictionary
    for atom in neighbors
]

neighbor_colors = [
    atom_colors.get(atom.element, "magenta")  # Use a default color if element is not in the dictionary
    for atom in neighbors
]

# Calculate distances from the atom of interest to each neighbor
distances = np.linalg.norm(neighbors.positions - atom_of_interest.position, axis=1)

# Calculate steric hindrance (number of overlaps based on van der Waals radii)
steric_score = sum(1 for d, r in zip(distances, neighbor_radii) if d < r)

# Output results
print(f"Atom of interest: {atom_of_interest}")
print(f"Number of neighbors within {cutoff_radius} Å: {len(neighbors)}")
print(f"Steric hindrance score: {steric_score}")

# Optional: Output detailed information about neighbors
for neighbor, distance, radius in zip(neighbors, distances, neighbor_radii):
    print(f"Neighbor: {neighbor}, Distance: {distance:.2f} Å, Radius: {radius:.2f} Å")

# --- Plotting Section ---
fig, ax = plt.subplots(figsize=(8, 8))

# Get positions of the atom of interest and neighbors
atom_pos = atom_of_interest.position
neighbor_positions = neighbors.positions

# Plot all atoms except the atom of interest (without legend)
for atom in u.atoms:
    if atom.index == atom_of_interest.index:
        continue
    color = atom_colors.get(atom.element, "gray")  # Default to gray if not defined
    ax.scatter(atom.position[0], atom.position[1], color=color, alpha=0.6, s=50)  # No label for other atoms

# Plot the atom of interest (in black)
ax.scatter(atom_pos[0], atom_pos[1], color='black', label='Atom of Interest (O)', s=100)

# Track labels to avoid duplicates in the legend
seen_labels = set()

# Plot neighbors and their circles
for neighbor, neighbor_pos, radius, color in zip(neighbors, neighbor_positions, neighbor_radii, neighbor_colors):
    label = f"Neighbor ({neighbor.element})"
    if label not in seen_labels:  # Add label only if not already added
        ax.scatter(neighbor_pos[0], neighbor_pos[1], color=color, label=label, s=80)
        seen_labels.add(label)
    else:
        ax.scatter(neighbor_pos[0], neighbor_pos[1], color=color, s=80)
    
    circle = plt.Circle(neighbor_pos[:2], radius, color=color, fill=False, linestyle='--', alpha=0.5)
    ax.add_artist(circle)

# Draw the cutoff circle
cutoff_circle = plt.Circle(atom_pos[:2], cutoff_radius, color='green', fill=False, linestyle='-', linewidth=2, label='Cutoff Radius')
ax.add_artist(cutoff_circle)

# Labeling and annotations
ax.set_title("2D Steric Hindrance Visualization")
ax.set_xlabel("X-coordinate (Å)")
ax.set_ylabel("Y-coordinate (Å)")
ax.legend(loc='upper right', fontsize=8)
ax.grid(True)
ax.set_aspect('equal', adjustable='box')

# Set limits based on the positions
padding = 1.0
x_min, x_max = atom_pos[0] - cutoff_radius - padding, atom_pos[0] + cutoff_radius + padding
y_min, y_max = atom_pos[1] - cutoff_radius - padding, atom_pos[1] + cutoff_radius + padding
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Show the plot
plt.show()

