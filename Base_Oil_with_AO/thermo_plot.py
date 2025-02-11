# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Feb 26 10:18:51 2024
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
plt.rc('font', size=18) # Default text sizes
plt.rc('axes', titlesize=12) # Axes title size
plt.rc('axes', labelsize=20) # Axes label size
plt.rc('xtick', labelsize=18) # X-axis tick label size
plt.rc('ytick', labelsize=18) # Y-axis tick label size
plt.rc('legend', fontsize=18) # Legend fontsize

directory = [r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_80_O2_Kowalik_2019\Equilibrate\200NVT+500NPT',
             r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_80_O2_6A_Kowalik_2019\Equilibrate\200NVT+500NPT']


labels = ['25 PAO+80 O$_2$', '25 PAO+6 A+80 O$_2$']
colors = ['tab:blue','tab:red','tab:grey']
for i, dirr in enumerate(directory):
    filename = '\\log.lammps'
    logfile = dirr+filename
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
    plt.plot(ps[start:end],thermo_property[start:end],color=colors[i],
             label=labels[i])
    
plt.ylabel(ne.getlab(plot_property))
plt.xlabel('Time (ps)')
plt.legend(loc='upper right')
# plt.title(directory[directory.find('Base'):]+'\n\n')
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\thermoplot\thermoplot-'+ne.randstr()+'.png', dpi=300,bbox_inches='tight')
print(ps[-1])