# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 28 04:58:56 2024
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import pandas as pd
import numpy as np

def generate_gaussian_opt_input(filename, atom_symbols, positions,
                        method="B3LYP",
                        basis_set="6-31G*",
                        charge=0,
                        multiplicity=1, nproc = 20, memory = "90Gb",
                        xyzfile=False):
    

    if os.path.exists(filename):
        print('A gaussian opt file with the same name already exists. Skipping!')
        return
    
    lines = [
        f"%NProcShared={nproc}",
        f"%mem={memory}",
        f"#n {method}/{basis_set} Opt freq",
        "",
        "optimization",
        "",
        f"{charge} {multiplicity}",
    ]
    
    positions = np.array(positions)
    for atom, position in zip(atom_symbols, positions):
        lines.append(f"{atom} {position[0]:.6f} {position[1]:.6f} {position[2]:.6f}")
    lines.append("\n")
    
    with open(filename, "w") as f:
        f.write("\n".join(lines))
    
    print("Gaussian opt input file generated.")
    
    if xyzfile:
        write_xyzfile(filename+'.xyz', atom_symbols, positions)


def extract_optimized_coordinates(input_file):
    # Flags and data storage
    optimization_found = False
    coordinates_start = False
    coordinates_data = []

    with open(input_file, 'r') as file:
        for line in file:
            # Check for the "Optimization completed" section
            if 'Optimization completed' in line:
                optimization_found = True

            # If optimization is found, look for the coordinates block
            if optimization_found and 'Standard orientation:' in line:
                # Reset the coordinates data for the new block
                coordinates_data = []
                # Skip the next 4 lines (headers and separators)
                for _ in range(4):
                    next(file)
                coordinates_start = True
                continue

            # Collect coordinates if in the section
            if coordinates_start:
                if '-----' in line:  # End of coordinates section
                    coordinates_start = False
                else:
                    parts = line.split()
                    if len(parts) >= 6:  # Ensure line has the correct format
                        center = int(parts[0])
                        atomic_number = int(parts[1])
                        x = float(parts[3])
                        y = float(parts[4])
                        z = float(parts[5])
                        coordinates_data.append([center, atomic_number, x, y, z])

    # Convert to pandas DataFrame
    df = pd.DataFrame(coordinates_data, columns=['Center', 'Atomic Number', 'X', 'Y', 'Z'])

    return df

def write_xyzfile(filename, atom_symbols, positions):
    if os.path.exists(filename):
        print('A xyz file with the same name already exists. Skipping!')
        return
    
    lines=[f"{len(atom_symbols)}\n"]
    positions = np.array(positions)
    for atom, position in zip(atom_symbols, positions):
        lines.append(f"{atom} {position[0]:.6f} {position[1]:.6f} {position[2]:.6f}")
    
    # Step 4: Write to a file
    with open(filename, "w") as f:
        f.write("\n".join(lines))
    
    print("XYZ file generated.")