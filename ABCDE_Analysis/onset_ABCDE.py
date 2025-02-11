# -*- coding: utf-8 -*-
"""
Created on Fri May 26 00:12:36 2023

@author: shihab
"""
import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
#%%
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\C\C_300_O2\Production\1277K'

timestep = 0.25
index = directory.find('ABCDE')
x = directory[index+6:index+7]
ramp_rate      = 4
initial_temp   = 300
print('Molecule:',x)
print('Ramp Rate',ramp_rate,'K/ps')
print('Initial Temperature:',initial_temp)

mean_temp = np.array([])
mean_sc   = np.array([])
onsets    = []
tempss    = []
scss      = []

for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    print('------------',sim[1:],'------------')
    d = directory+sim
    filename = '\\bonds.reaxc'
    bfpath   = d+filename # bond file path
    if not os.path.exists(bfpath):
        continue
    atominfo = bfp.get_neighbours(bfpath)
    neighbours,atomtypes = atominfo
    atomsymbols = ['H','C','O']
    lignins = {'A':'H12C10O3',
               'B':'H34C26O4',
               'C':'H44C29O2',
               'D':'H30C19O3',
               'E':'H14C11O4',
               'e':'H62C30'}
    
    main_molecule = lignins[x]
    
    temp,sc = bfp.get_speciesVStemp(main_molecule, *atominfo, atomsymbols,it=initial_temp,ramp=ramp_rate,ts=timestep)
    onset,_   = bfp.get_onset(temp, sc,imc=50,ig=[1700,107,0,50])
    tempss.append(temp)
    scss.append(sc)
    onsets.append(onset)
    print('Onset:',onset,'\n')
#%%
mean_temp         = sum(tempss)/len(tempss)
mean_sc           = sum(scss)/len(scss)

mean_onset,sc_fit = bfp.get_onset(mean_temp, mean_sc, imc=50,ig=[1100,107,0,50])
print('Onsets: {:.0f}, {:.0f}, {:.0f}'.format(*onsets))
print('Mean Onset: {:.0f}'.format(mean_onset))
#finding error
error = 0
for onset in onsets:
    error += (onset-mean_onset)**2
error = (error/len(onsets))**0.5
print('error: {:.1f}'.format(error))

#--------plotting-----
d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n\nMolecule {}\n'.format(x)+'Onset {:.0f} $\pm$ {:.1f} K\n'.format(mean_onset,error)

plt.style.use('classic')
plt.scatter(mean_temp,mean_sc,marker='s', label='data',color='#77dd77',s=50,alpha=1)
plt.plot(mean_temp,sc_fit, label='fit',color='r',linewidth=2)
plt.plot([mean_onset]*1000,np.linspace(-5,25,1000),'--',color='black',linewidth=2)
plt.title(title)
plt.ylim(-1,26)
plt.xlim(00,3550)
plt.xlabel('Temperature (K)')
plt.ylabel('Number of molecule')
plt.legend()
#plt.grid(visible=True)
plt.savefig('..\\python_outputs\\onset\\onset-'+str(random.randint(1000000,9999999)), dpi=400,bbox_inches='tight')
plt.show()
#-------------------