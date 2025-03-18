# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Mar 11 01:40:16 2025
"""

import numpy as np

# Predefined van der Waals radii (in Angstroms)
VDW_RADII = {
    1: 1.70,   # Carbon (C) - Mass 12.01115
    2: 1.70,   # Carbon (C) - Mass 12.01115
    3: 1.70,   # Carbon (C) - Mass 12.01115
    4: 1.70,   # Carbon (C) - Mass 12.01115
    5: 1.70,   # Carbon (C) - Mass 12.01115
    6: 1.75,   # Chlorine (Cl) - Mass 35.453
    7: 1.70,   # Carbon (C) - Mass 12.01115
    8: 1.20,   # Hydrogen (H) - Mass 1.00797
    9: 1.20,   # Hydrogen (H) - Mass 1.008
    10: 1.55,  # Nitrogen (N) - Mass 14.0067
    11: 1.52,  # Oxygen (O) - Mass 15.9994
    12: 1.52,  # Oxygen (O) - Mass 15.9994
    13: 1.52,  # Oxygen (O) - Mass 15.9994
    14: 1.52   # Oxygen (O) - Mass 15.9994
}

def read_lammps_data(lammps_file):
    """
    Reads atomic information from a LAMMPS data file.

    Parameters:
        lammps_file (str): Path to the LAMMPS data file.

    Returns:
        list of tuples: (atom_type, x, y, z)
    """
    with open(lammps_file, 'r') as file:
        lines = file.readlines()

    atom_section_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atom_section_start = i + 2  # Skip header
            break

    if atom_section_start is None:
        raise ValueError("No 'Atoms' section found in the LAMMPS data file.")

    atom_data = []
    for line in lines[atom_section_start:]:
        if line.strip() == "":
            break
        parts = line.split()
        atom_type = int(parts[2])  # Atom type (LAMMPS format)
        x, y, z = map(float, parts[4:7])  # Coordinates
        atom_data.append((atom_type, x, y, z))

    return atom_data

def compute_occupied_volume(lammps_file):
    """
    Calculates the occupied volume using van der Waals radii.

    Parameters:
        lammps_file (str): Path to the LAMMPS data file.

    Returns:
        float: Total occupied volume in cubic Angstroms.
    """
    atoms = read_lammps_data(lammps_file)

    occupied_volume = 0.0
    for atom_type, _, _, _ in atoms:
        vdw_radius = VDW_RADII.get(atom_type, 1.50)  # Default to 1.50 Å if not found
        print(atom_type, vdw_radius)
        atom_volume = (4/3) * np.pi * (vdw_radius ** 3)
        occupied_volume += atom_volume

    return occupied_volume

# Example Usage
lammps_file = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\Chloroethane\PCFF-Gastiger\Eq\Bulk_Phase\Sim-1\eq.npt-density.data"  # Replace with your actual file path
V_occ = compute_occupied_volume(lammps_file)
print(f"Estimated Occupied Volume: {V_occ:.2f} Å³")
