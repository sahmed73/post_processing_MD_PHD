# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 29 00:12:42 2024
"""

import cclib
from rdkit.Chem import GetPeriodicTable
import os

def gaussian_to_xyz(filename, output_xyz):
    """
    Convert Gaussian optimization output to XYZ trajectory file.

    Parameters:
    - filename: Path to the Gaussian .log or .out file.
    - output_xyz: Path to save the resulting XYZ file.
    """
    # Parse Gaussian output using cclib
    data = cclib.io.ccread(filename)

    # Extract atom numbers and optimization coordinates
    atom_types = data.atomnos  # Atomic numbers
    coordinates = data.atomcoords  # All geometry steps

    # Get periodic table for atomic symbols
    pt = GetPeriodicTable()

    # Write to XYZ file
    with open(output_xyz, "w") as f:
        for step, coord in enumerate(coordinates):
            f.write(f"{len(atom_types)}\n")
            f.write(f"Step {step + 1}\n")
            for atom, xyz in zip(atom_types, coord):
                symbol = pt.GetElementSymbol(int(atom))  # Convert atomic number to element symbol
                f.write(f"{symbol} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")
    print(f"XYZ trajectory saved to: {output_xyz}")


# Example usage
dirr=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\A0003\pre'
filename = dirr+"\\opt.out"  # Replace with your Gaussian .log or .out file
output_xyz = dirr+"\\trajectory.xyz"  # Replace with your desired output file name

gaussian_to_xyz(filename, output_xyz)
