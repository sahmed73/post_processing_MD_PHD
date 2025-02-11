# -*- coding: utf-8 -*-
"""
Created on Sun May 21 21:24:29 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
#%%
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_100_O2\Onset'
timestep = 0.25

mean_temp = np.array([])
mean_sc   = np.array([])
onsets    = []
tempss    = []
scss      = []

for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    d = directory+sim
    filename = '\\bonds.reaxc'
    bfpath   = d+filename # bond file path
    if not os.path.exists(bfpath):
        continue
    atominfo = bfp.get_neighbours(bfpath)
    neighbours,atomtypes = atominfo
    atomsymbols = ['H','C','O']
    main_molecule = 'H12C10O3'
    
    temp,sc = bfp.get_speciesVStemp(main_molecule, *atominfo, atomsymbols,it=300,ramp=4.0)
    onset,_   = bfp.get_onset(temp, sc,imc=50,show_fit=True)
    tempss.append(temp)
    scss.append(sc)
    onsets.append(onset)
#%%
mean_temp         = sum(tempss)/len(tempss)
mean_sc           = sum(scss)/len(scss)

mean_onset,sc_fit = bfp.get_onset(mean_temp, mean_sc, imc=50)
print(onsets,mean_onset)
#finding error
error = 0
for onset in onsets:
    error += (onset-mean_onset)**2
error = (error/len(onsets))**0.5
print('error:',error)

#--------plotting-----
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\'):]+'\n'+'Onset {:.0f} $\pm$ {:.1f} K\n'.format(mean_onset,error)

plt.style.use('classic')
plt.scatter(mean_temp,mean_sc,marker='s', label='data',color='dimgrey',s=50,alpha=1)
plt.plot(mean_temp,sc_fit, label='fit',color='r',linewidth=2)
plt.plot([mean_onset]*1000,np.linspace(-5,50,1000),'--',color='black',linewidth=2)
plt.title(title)
#plt.errorbar(ex, ey, yerr=estd,fmt='o',color='black')
#plt.xticks(np.arange(1000, 2600, 500))
plt.ylim(-1,53)
plt.xlim(0,3550)
plt.xlabel('Temperature (K)')
plt.ylabel('Number of molecule')
plt.legend()
#plt.grid(visible=True)
plt.savefig('python_outputs\\onset\\onset-'+str(random.randint(999999, 1000000)), dpi=400,bbox_inches='tight')
plt.show()
#-------------------