# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:24:40 2023
Last checked on 12/10/2023

@author: sahmed73

## Purpose:
This script reads a species output file generated by the `reaxff/species` command in LAMMPS 
and converts it into a Pandas DataFrame with time on one axis and species counts on the other.

## Expected Species File Format:
The input file should have repeating blocks where each block starts with a header line containing:
- `Timestep`, `No_Moles`, `No_Specs`, followed by species names (e.g., `H62C30`, `O2`)
- The next line contains numerical values for the corresponding timestep.

Example format:

# Timestep     No_Moles     No_Specs     H62C30	O2	
1000        225          2	 25	 200	
# Timestep     No_Moles     No_Specs     H62C30	O2	
2000        225          2	 25	 200	
# Timestep     No_Moles     No_Specs     H62C30	O2	
3000        225          2	 25	 200	
# Timestep     No_Moles     No_Specs     H62C30	O2	
4000        225          2	 25	 200	

"""

import pandas as pd

def get_species(speciesfile, timestep=None):
    """
    Reads a species output file and converts it into a Pandas DataFrame.

    Parameters:
    - speciesfile (str): Path to the species output file.
    - timestep (float, optional): The time step size in femtoseconds (fs).
      If provided, the function converts timesteps to picoseconds (ps).

    Returns:
    - df (pd.DataFrame): A DataFrame where:
      - Rows represent timesteps (or time in ps if timestep is given).
      - Columns represent species counts.

    Usage Example:
    df = get_species("species.out", timestep=0.25)
    """

    with open(speciesfile, 'r') as sf:
        species = {}  # Dictionary to store species counts for each timestep
        for line in sf:
            if "Timestep" in line:  # Identify the header line
                headers = line.strip().split()[1:]  # Ignore the '#' symbol
                species_name = headers[3:]  # Extract species names (skip first 3 columns)
            else:
                values = [int(x) for x in line.strip().split()]  # Convert values to integers
                step = values[0]  # Extract timestep value
                species_count = values[3:]  # Extract species counts (after No_Moles, No_Specs)
                
                species[step] = {}
                for key, value in zip(species_name, species_count):
                    if key in species[step]:
                        species[step][key] += value  # Accumulate count if species already exists
                    else:
                        species[step][key] = value  # Initialize species count
    
    # Convert dictionary to DataFrame
    df = pd.DataFrame(species).fillna(0).T  # Fill NaN with 0 and transpose for correct format
    df.index.name = 'Timestep'

    # Convert timesteps to time (in picoseconds) if `timestep` is provided
    if timestep is not None:
        df['Time'] = df.index * timestep / 1000  # Convert fs to ps
        df.index = df['Time']  # Replace index with time in ps
        df.drop(['Time'], axis=1, inplace=True)  # Remove extra column

    return df


#%% User Input Section
dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"

# Path to the species output file
speciesfile_path = dirr + r"\\species.out"

df = get_species(speciesfile_path, timestep=0.25)  # Here, timestep is 0.25 fs

print(df)


#%% Plotting Section
import matplotlib.pyplot as plt


plt.plot(df['H62C30'], label='$H_{62}C_{30}$', color='k')
plt.plot(df['O2'], label='$O_{2}$', color='red')
plt.xlabel('Time (ps)')
plt.ylabel('Number of Species')
plt.legend()
plt.show()
