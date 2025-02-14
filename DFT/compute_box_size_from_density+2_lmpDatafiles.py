# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Feb 13 12:16:00 2025
"""

import numpy as np
from periodictable import elements

def read_lammps_data(filename):
    """Extract atomic symbols and molecule composition from a LAMMPS data file."""
    atom_symbols = []
    with open(filename, 'r') as file:
        inside_atoms_section = False
        for line in file:
            if line.strip().startswith("Atoms"):
                inside_atoms_section = True
                continue
            if inside_atoms_section and line.strip() == "":
                continue # skip the blank lines
            
            if inside_atoms_section and line[0].isalpha():
                break  # Stop when leaving Atoms section

            if inside_atoms_section:
                parts = line.split("#")
                if len(parts) > 1:
                    symbol = parts[1].strip().split('/')[-1]  # Extract element after "/"
                    atom_symbols.append(symbol)

    return atom_symbols

def molecular_weight(symbols):
    """Calculate molecular weight from atomic symbols."""
    return sum(elements.symbol(sym).mass for sym in symbols)

def compute_box_size(lammps_files, molecule_counts, target_density):
    """
    Compute cubic box size given multiple LAMMPS data files, molecule counts, and density.

    Parameters:
    - lammps_files: List of LAMMPS data filenames.
    - molecule_counts: List of molecule counts (same order as files).
    - target_density: Target density in g/cm³.

    Returns:
    - Box size in Å.
    """
    total_mass = 0

    for i, lammps_file in enumerate(lammps_files):
        symbols = read_lammps_data(lammps_file)
        mol_weight = molecular_weight(symbols)
        print(mol_weight)
        total_mass += mol_weight * molecule_counts[i]

    # Convert density from g/cm³ to atomic mass unit (amu) per Å³
    density_in_amu_per_a3 = target_density * 0.6022140857  # Conversion factor

    # Compute volume in Å³
    volume = total_mass / density_in_amu_per_a3

    # Compute cubic box length
    box_length = volume ** (1/3)

    return box_length

# Example usage
lammps_files = [r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0002\DataFile\from_LUNAR\from_bond_react_merge\A0002_pre_typed_IFF_merged.data",
                r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0002\DataFile\from_LUNAR\from_bond_react_merge\PAOr_pre_typed_IFF_merged.data"]  # Replace with actual filenames
molecule_counts = [50, 100]  # Adjust based on actual molecule numbers
target_density = 0.20  # g/cm³

box_size = compute_box_size(lammps_files, molecule_counts, target_density)
print(f"Computed cubic box size: {box_size:.3f} Å")
