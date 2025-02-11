# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 06:16:59 2023

@author: arup2
"""

import magnolia.log_parser_FHB as lfp
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import magnolia.needless_essential as ne
import math
import random

from scipy.optimize import curve_fit
plt.style.use('classic')

base_oil = ["PAO4","Squalane"]
color    = ["blue","green"]

for i,oil in enumerate(base_oil):
    directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\Only_25_PAO4_New\Equilibrate\Different_Potentials\200NVT+1000NPT_CHO_2016'
    
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
        
    ps = ps - 200
    
    thermo_property = np.array(thermo[plot_property])
    if plot_property in ['PotEng','TotEng','KinEng']:
        thermo_property = thermo_property - min(thermo_property)
    
    if plot_property in ['Volume']:
        thermo_property /= 1000
        ylabel = 'Volume (nm$^3$)'
    else:
        ylabel = ne.getlab(plot_property)
     
    start = bfp.get_nearestindex(ps, key='first')
    end = bfp.get_nearestindex(ps, key='last')
    
    plt.plot(ps[start:end],thermo_property[start:end],color='r',label=oil)
    
    if showstep:
        plt.xlabel('Iteration')
    else:
        plt.xlabel(ne.getlab('time'))
    
    
    plt.ylabel(ylabel)
    plt.title(directory[directory.find('Base_Oil'):]+'\n\n')
    plt.ylim(0.6,1.2)
    plt.xlim(0,400)
    # plt.legend()
    # plt.plot([700]*10,np.linspace(0,1.2,10),'--')
    plt.savefig('..\\python_outputs\\thermoplot\\thermoplot_{}'.format(random.randint(0,1000000000000)), dpi=300,bbox_inches='tight')


