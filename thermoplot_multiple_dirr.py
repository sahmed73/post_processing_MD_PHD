# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Mar  1 05:35:32 2024
"""


import magnolia.log_parser_FHB as lfp
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import magnolia.needless_essential as ne
import math

from scipy.optimize import curve_fit
plt.style.use('classic')
# Set default font sizes
plt.rc('font', size=10) # Default text sizes
plt.rc('axes', titlesize=10) # Axes title size
plt.rc('axes', labelsize=20) # Axes label size
plt.rc('xtick', labelsize=15) # X-axis tick label size
plt.rc('ytick', labelsize=15) # Y-axis tick label size
plt.rc('legend', fontsize=15) # Legend fontsize


common = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4'

directory = [r'\25_PAO4_80_O2_Soria\Equilibrate\200NVT+500NPT',
             r'\25_PAO4_80_O2_6A_Soria\Equilibrate\200NVT+500NPT',
             r'\20_PAO4_80_O2_6A_Soria\Equilibrate\200NVT+500NPT']

filename = '\\log.lammps'

labels = ['25 PAO, 80 O$_2$',
          '25 PAO, 80 O$_2$, 6 A',
          '20 PAO, 80 O$_2$, 6 A']

for i, dirr in enumerate(directory):
    logfile = common+dirr+filename
    timestep = 0.25
    thermo = lfp.thermo_dict_v2(logfile,serial=2)
    
    plot_property = 'PotEng'
    showstep = False
    steps = thermo['Step']
    
    
    if showstep:
        ps = np.array(steps)
    else:
        ps = (np.array(bfp.step2picosecond(steps, timestep)))
        
    
    thermo_property = np.array(thermo[plot_property])
    if plot_property in ['PotEng','TotEng','KinEng']:
        thermo_property = thermo_property - min(thermo_property)
     
    start = 0
    end = -1
    
    ps = ps-ps[0]
    plt.plot(ps[start:end],thermo_property[start:end],
             label=labels[i])
    
    if showstep:
        plt.xlabel('Iteration')
    else:
        plt.xlabel(ne.getlab('time'))
        
plt.ylabel(ne.getlab(plot_property))
# plt.title(labels)
plt.legend(loc='upper right')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\thermoplot\thermoplot-'+ne.randstr()+'.png', dpi=300,bbox_inches='tight')
print(ps[-1])