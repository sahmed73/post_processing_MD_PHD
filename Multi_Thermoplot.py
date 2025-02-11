# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 08:30:48 2023

@author: arup2
"""

import magnolia.log_parser_FHB as lfp
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import magnolia.needless_essential as ne
import math

from scipy.optimize import curve_fit
plt.style.use('classic')

dirs = {'PAO4':r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\Equilibrate\200NVT_500NPT_500NPT_Combined\No_Velocity',
        'Squalane':r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\Equilibrate\200NVT_500NPT_500NPT_Combined\No_Velocity'}

for baseoil in ['PAO4','Squalane']:
    directory = dirs[baseoil]
    filename = '\\log.lammps'
    logfile = directory+filename
    timestep = 0.25
    thermo = lfp.thermo_dict(logfile,3)
    
    
    plot_property = 'Density'
    showstep = False
    
    steps = thermo['Step']
    if showstep:
        ps = np.array(steps)
    else:
        ps = (np.array(bfp.step2picosecond(steps, timestep)))
    
    ps = ps-200
    
    thermo_property = np.array(thermo[plot_property])
    if plot_property in ['PotEng','TotEng','KinEng']:
        thermo_property = thermo_property - min(thermo_property)
     
    start = bfp.get_nearestindex(ps, key='first')
    end = bfp.get_nearestindex(ps, key='last')
    
    plt.plot(ps[start:end],thermo_property[start:end],label=baseoil)

if showstep:
    plt.xlabel('Iteration')
else:
    plt.xlabel(ne.getlab('time'))
    
plt.ylabel(ne.getlab(plot_property))
plt.title(directory[directory.find('Base_Oil'):]+'\n\n')
plt.xlim(0,1005)
plt.legend(loc='upper left')
plt.savefig('python_outputs\\thermoplot\\thermoplot-'+ne.randstr(), dpi=300,bbox_inches='tight')

# need iteration for a temperature
d = dict(zip(thermo_property,ps))
onset = 1420
for key,value in d.items():
    if value%4000==0 and abs(key-onset)<10:
        print('({},{})'.format(key,value))
# finding restart file

