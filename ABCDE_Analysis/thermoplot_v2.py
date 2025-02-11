# -*- coding: utf-8 -*-
"""
Created on Sun May 28 01:45:11 2023

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

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\200_Less_TempRamp_1150K'
filename = '\\log.lammps'
logfile = directory+filename
timestep = 0.25
thermo = lfp.thermo_dict_v2(logfile,serial=2)

plot_property = 'Temp'
showstep = False

steps = thermo['Step']
if showstep:
    ps = np.array(steps)
else:
    ps = (np.array(bfp.step2picosecond(steps, timestep)))
    

thermo_property = np.array(thermo[plot_property])
if plot_property in ['PotEng','TotEng','KinEng']:
    thermo_property = thermo_property - min(thermo_property)
 
start = bfp.get_nearestindex(ps, key='first')
end = bfp.get_nearestindex(ps, key='last')

plt.plot(ps[start:end],thermo_property[start:end],color='black')

if showstep:
    plt.xlabel('Iteration')
else:
    plt.xlabel(ne.getlab('time'))
    
#plt.yaxis.set_ticks(np.arange(0, 85, step=20))
# plt.text(0,800,'<-------NPT------->',fontsize=18,color='red')
# plt.text(300,800,'|<--------------NVT-------------->',fontsize=18,color='red')
#plt.ylim(0.5,1.1)
plt.ylabel(ne.getlab(plot_property))
plt.title(directory[directory.find('Antioxidants'):]+'\n\n')
plt.savefig('..\\python_outputs\\thermoplot\\thermoplot-'+ne.randstr(), dpi=300,bbox_inches='tight')