# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Mar  6 20:26:14 2024
"""

import pandas as pd

def num_and_types(filename, data_lines):
    # Read the specified lines from the file. Note that Pandas uses zero-based indexing.
    # MATLAB's data_lines [start, end] is inclusive, so adjust the end index for Python.
    data = pd.read_csv(filename, 
                       skiprows=data_lines[0]-1, 
                       nrows=data_lines[1]-data_lines[0]+1, 
                       delim_whitespace=True, 
                       header=None, 
                       usecols=[0],  # Assuming 'LAMMPS' column corresponds to the first column
                       names=['LAMMPS'])
    
    # Convert the DataFrame to a numpy array and return
    return data['LAMMPS'].values

def import_dims(filename, data_lines):
    # Read the specified lines from the file
    # Adjust for zero-based indexing in Python and read the specific lines and columns
    data = pd.read_csv(filename, 
                       skiprows=data_lines[0]-1, 
                       nrows=data_lines[1]-data_lines[0]+1, 
                       delim_whitespace=True, 
                       header=None, 
                       usecols=[0, 1],  # Assuming 'LAMMPS' is the first column and 'data' is the second
                       names=['LAMMPS', 'data'])
    
    # Convert the selected columns to a numpy array and return
    return data.values

def import_masses(filename, data_lines):
    # Read the specified lines from the file for the 'masses' section
    data = pd.read_csv(filename, 
                       skiprows=data_lines[0]-1, 
                       nrows=data_lines[1]-data_lines[0]+1, 
                       delim_whitespace=True, 
                       header=None, 
                       usecols=[0, 1],  # Assuming 'LAMMPS' and 'data' are the first two columns
                       names=['LAMMPS', 'data'])
    
    # Convert the DataFrame to a numpy array and return
    return data.values

def import_pair_coeffs(filename, data_lines):
    # Read the specified lines from the file for the 'pairCoeffs' section
    data = pd.read_csv(filename, 
                       skiprows=data_lines[0]-1, 
                       nrows=data_lines[1]-data_lines[0]+1, 
                       delim_whitespace=True, 
                       header=None, 
                       usecols=[0, 1, 2],  # Adjust based on actual column positions
                       names=['LAMMPS', 'data', 'file'])
    
    # Convert the DataFrame to a numpy array and return
    return data.values

# Example usage
filename = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\OPLSAA\PAO+Molecule-A\PAO4.lmp'
data_lines = [3, 7]  # Example line numbers for the section to import
num = num_and_types(filename, data_lines)
print(num)
print('----------')

data_lines = [15, 17]  # Example line numbers for the section to import
dims = import_dims(filename, data_lines)
print(dims)
print('----------')

data_lines = [3, 7]  # Adjust these line numbers based on your file
masses = import_masses(filename, data_lines)
print(masses)
print('----------')