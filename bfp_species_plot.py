# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 02:40:31 2023

@author: Shihab
"""

'''
bond file formate
id type nb id_1 .... id_nb mol bo_1 ... bo_nb abo nlp q
'''
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import re
import magnolia.bondfile_parser as bfp
import random

#%%---------------------main-----------------------------------
start_time = time.time()
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_100_O2\Production\200K_less'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
neighbours,atomtype = bfp.get_neighbours(bondfile,atypes=True)


#%%----------------------------------------------------------
AO_1 = 'H12C10O3' 
AO_2 = 'H34C26O4'
seek_species =  ['H2CO','H4CO','H2O','CO2']#
label = bfp.make_molecular_formula_latex(seek_species,sort=True) 
color = ['red','green','blue','black']
steps = list(neighbours.keys())
atomsymbols = 'HCO'

x_list = []
y_list = []
for j in range(len(seek_species)):
    nspecies = [0]*len(steps)
    for i in range(len(steps)):
        molecules = bfp.get_molecules(neighbours[steps[i]])
        speciesList = bfp.get_species(molecules, atomtype, atomsymbols)
        species = seek_species[j]
        if species in speciesList.keys():
            nspecies[i] = speciesList[species]

    x_list.append(steps.copy())
    y_list.append(nspecies.copy())
#%%----------------------------------------------------------

for i in range(len(seek_species)):
    x = bfp.step2picosecond(x_list[i], 0.25)
    y = y_list[i]
    start = bfp.get_nearestindex(x, 100)
    end   = -1#get_timeindex(x, 110)
    plt.plot(x[start:end],y[start:end],label=label[i],color=color[i],linewidth=1.0)
    
plt.xlabel('Time (ps)')
plt.ylabel('Number of molecules')
plt.yticks(np.arange(0, 84, step=5))
#plt.ylim(41)
plt.legend()
plt.title(directory[directory.find('AO_Oxidation'):]+'\n\n')
rand = str(random.randint(100000, 999999))
plt.savefig('python_outputs//figures//species-time-plot-'+rand,dpi=300,bbox_inches='tight')

runtime = time.time()-start_time
print('----Total run Time: {} min {:.0f} sec---------'.format(int(runtime/60),runtime%60))
################### Beep after finish running #########################
'''winsound.Beep(1200, 100)
winsound.Beep(2000, 100)
winsound.Beep(1200, 100)'''