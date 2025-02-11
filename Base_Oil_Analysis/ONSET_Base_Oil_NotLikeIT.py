# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 12 23:13:41 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd
from scipy import optimize


### Directory ###
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Onset"
filename  = "\\bonds.reaxc"


whole = 'H62C30'
whole_count = {}
for sim in ["Sim-1","Sim-2","Sim-3"]:
    location = directory+"\\"+sim
    bondfilepath = location+filename
    bonddata = bfp.parsebondfile(bondfilepath)
    
    neighbours = bonddata['neighbours']
    atypes     = bonddata['atypes']
    asyms      = ['H','C','O']
    steps      = np.array(list(neighbours.keys()))
    
    whole_count[sim] = []
    for step, neigh in neighbours.items():
        molecules = bfp.get_molecules(neigh)
        count = 0
        for molecule in molecules:
            species = bfp.get_molecular_formula(molecule, atypes, asyms)
            if species == whole:
                count+=1
        whole_count[sim].append(count)

#%%
df = pd.DataFrame(whole_count)
print(df)


# Melting the DataFrame to Long Format
df_melted = pd.melt(df.reset_index(), id_vars=['index'], 
                    value_vars=['Sim-1', 'Sim-2', 'Sim-3'],
                    var_name='Simulation', value_name='Molecule Count')

# Curve Fitting Example: (Assuming a linear fit for simplicity)
def model_func(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_fit(y,A,B,C,D):
    x = A + B*np.log((D-y)/(y-C))
    return x

# Concatenating all data (or average them, depending on your fit approach)
all_y = pd.concat([df[col] for col in df.columns])
all_x = np.tile(df.index, df.shape[1])
popt, _ = optimize.curve_fit(model_func, all_x, all_y)

# Generating Data for Fitted Curve
x_line = np.arange(min(df.index), max(df.index), 1)
y_line = model_func(x_line, *popt)

# Plotting Data and Fitted Curve
# sns.scatterplot(data=df_melted, x='index', y='Molecule Count', hue='Simulation', palette='viridis')
plt.plot(x_line, y_line, '--', color='red', label='Fitted Curve')

# Shading (using a computed upper and lower bound or min-max values)
lower_bound = df.min(axis=1)
upper_bound = df.max(axis=1)
plt.fill_between(df.index, lower_bound, upper_bound, color='gray', alpha=0.4)

plt.legend()
plt.title('Molecule Count Over Time with Fitted Curve')
plt.xlabel('Time')
plt.ylabel('Molecule Count')
plt.show()

