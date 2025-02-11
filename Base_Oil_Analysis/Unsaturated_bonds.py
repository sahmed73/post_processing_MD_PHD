# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 13 16:49:29 2023
"""
import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import magnolia.plot_template as mplt
from scipy.optimize import curve_fit
from scipy.optimize import fsolve

fig, ax = plt.subplots()
double = 1.2 # double bond
timestep = 0.25
initial_temp  = 300
temp_ramprate = 4
colors = ['gold', 'green']
oils = ['PAO4','Squalane']
simulations = ["Sim-1","Sim-2", "Sim-3"]

data = []
for oil in oils:
    df = pd.DataFrame()
    for sim in simulations:
        print(oil,sim)
    ### Directory ###
        directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(oil,oil,sim)
        filename  = "\\bonds.reaxc"
        bondfilepath = directory+filename
        
        ### Parsing Bondfile ###
        bonddata = bfp.parsebondfile(bondfilepath,bo=True)
        
        ## Geting neighbour lists
        neighbours = bonddata['neighbours']
        atypes     = bonddata['atypes']
        bondorders = bonddata['bondorders']
        asyms      = ['H','C','O']
        steps      = np.array(list(neighbours.keys()))
        
        ### Loop through the neighbours ###
        COBond_count = []
        for step, neigh in neighbours.items():
            count = 0
            for parent, children in neigh.items():
                for child in children:
                    bo = bondorders[step][parent][child]
                    if bo>double and atypes[parent]==2 and atypes[child]==2:
                        count+=1
            COBond_count.append(count)
        
        df[sim]=COBond_count
        
    df.index=(steps*timestep/1000)*temp_ramprate+initial_temp
    data.append(df)

#%% ploting using own function
# Composite model: logistic plus polynomial
# Set default font sizes
plt.rcParams['figure.figsize'] = [6, 4]
plt.rc('font', size=14) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=17) # Axes label size
plt.rc('xtick', labelsize=16) # X-axis tick label size
plt.rc('ytick', labelsize=16) # Y-axis tick label size
plt.rc('legend', fontsize=16) # Legend fontsize

colors = ['gold', 'green']
oils = ['PAO4','Squalane']
def composite_model(x, L, k, x0, a, b, c):
    # Logistic part
    logistic_part = L / (1 + np.exp(-k * (x - x0)))
    # Polynomial part (second degree)
    poly_part = a * x**2 + b * x + c
    # Combine the two parts
    return logistic_part + poly_part

def find_x_from_y(y_target, popt, x_guess):
    # Define the equation that fsolve will find the root for
    def equation(x):
        L, k, x0, a, b, c = popt
        logistic_part = L / (1 + np.exp(-k * (x - x0)))
        poly_part = a * x**2 + b * x + c
        return logistic_part + poly_part - y_target

    # Use fsolve to find the x value that gives the desired y value
    x_solution = fsolve(equation, x_guess)
    return x_solution[0]

fig, ax = plt.subplots(1,1)

for i, df in enumerate(data):
    adf = df.copy()
    adf.index = (df.index-1600)/4
    adf = adf.loc[0:,:]
    mplt.mean_range_plot(adf, color=colors[i],
                         alpha=0.1,linewidth=0.2, ax=ax)
    ax.plot([],[],color=colors[i],linewidth=3.5,label=oils[i])
        
    # x_data = df.index
    # y_data = df.mean(axis=1)
    # # Initial guesses for the composite model parameters
    # # For logistic part: L, k, x0
    # # For polynomial part: a, b, c
    # initial_guesses_composite = [np.max(y_data), 1, np.median(x_data),
    #                              0, 0, np.min(y_data)]
    
    # # Perform the curve fitting for the composite model
    # popt, _ = curve_fit(composite_model, x_data, y_data,
    #                                            p0=initial_guesses_composite,
    #                                            maxfev=8000)
        
    # # Plot the fitted curve
    # y_fit = composite_model(x_data, *popt)
    # # plt.plot(x_data, y_fit, color=colors[i], label='fit',linewidth=0.5)
    # onset = find_x_from_y(401, popt, x_guess=2000)
    # # plt.plot([onset]*2,[400,620],color=colors[i])
    # print(onset)

plt.legend(loc='upper left')
plt.xlabel('Time (ps)')
plt.ylabel('Number of unsatuareted bonds')
title = ",".join(oils+simulations)+'\n\n\n'
plt.title(title+'\n\n', fontsize=8)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\figures'
plt.savefig(savedir+'\\Unsaturated_bonds_{}'.format(random.randint(0, 10000000)), dpi=500, bbox_inches='tight')