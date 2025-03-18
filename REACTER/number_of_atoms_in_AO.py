# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  6 15:26:03 2024
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np

from ase import io
from ase.data import atomic_masses


def count_atoms_from_xyz(file_path):
    atoms = io.read(file_path, format='xyz')
    return len(atoms)

def calculate_molecular_weight(file_path):
    atoms = io.read(file_path, format='xyz')
    molecular_weight = sum(atomic_masses[atom.number] for atom in atoms)
    return molecular_weight

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size'] = 18

parent_dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001"
extensions = [r"\A0001\DataFile\Single_Molecules\A0001_pre.xyz",
              r"\A0002\DataFile\Single_Molecules\A0002_pre.xyz",
              r"\A0003\DataFile\Single_Molecules\A0003_pre.xyz",
              r"\A0004\DataFile\Single_Molecules\A0004_optimized_pre1.xyz",
              r"\A0005\DataFile\Single_Molecules\A0005_optimized_pre1.xyz"]

AOs = ['A0001', 'A0002', 'A0003', 'A0004', 'A0005']

for i, extension in enumerate(extensions):
    dirr = parent_dirr + extension

    xyz_file = dirr  # Replace with the actual file path
    num_atoms = count_atoms_from_xyz(xyz_file)
    mol_weight = calculate_molecular_weight(xyz_file)

    print(f"Number of atoms in {AOs[i]}: {num_atoms}")
    print(f"Molecular weight of {AOs[i]}: {mol_weight:.2f} g/mol")
