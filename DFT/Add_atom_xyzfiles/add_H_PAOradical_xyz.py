# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 29 01:30:46 2024
"""
# nooooooooooooooooooooooooooooooooooooooooooooooo

from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT  # Basic force-field calculator

def add_hydrogen_and_optimize(xyz_file, output_file, radical_index):
    """
    Add a hydrogen atom to a radical and optimize its geometry.

    Parameters:
    - xyz_file: Path to the input XYZ file.
    - output_file: Path to save the optimized XYZ file.
    - radical_index: Index (0-based) of the radical atom to which H will be added.
    """
    # Read the molecule
    from ase.io import read, write
    molecule = read(xyz_file)

    # Get the radical atom position
    radical_pos = molecule.positions[radical_index]

    # Estimate initial H position (bond length ~1 Å from the radical atom)
    h_pos = radical_pos + [1.0, 0.0, 0.0]  # Place H atom 1 Å away along x-axis

    # Add the hydrogen atom
    molecule.append('H')  # Add H to the molecule
    molecule.positions[-1] = h_pos  # Set H position

    # Set up a calculator for optimization
    molecule.set_calculator(EMT())  # Use EMT (simple force field)

    # Optimize geometry
    dyn = BFGS(molecule)
    dyn.run(fmax=0.01)  # Run until forces converge to below 0.01 eV/Å

    # Save the optimized structure
    write(output_file, molecule)
    print(f"Optimized structure with H saved to: {output_file}")

# Example usage
xyz_file = "radical.xyz"  # Input radical structure
output_file = "optimized_with_H.xyz"  # Output file with added H
radical_index = 0  # Index of the atom with the radical (adjust as needed)

add_hydrogen_and_optimize(xyz_file, output_file, radical_index)
