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

from scipy.optimize import curve_fit
plt.style.use('classic')

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600_0.1fs\Sim-1'

filename = '\\log.lammps'
logfile = directory+filename
timestep = 0.10
thermo = lfp.thermo_dict(logfile,2)

for pp in ['Temp']:
    plot_property = pp
    showstep = False
    
    steps = thermo['Step']
    if showstep:
        ps = np.array(steps)
    else:
        ps = (np.array(bfp.step2picosecond(steps, timestep)))
        
    ps = ps
    
    thermo_property = np.array(thermo[plot_property])
    if plot_property in ['PotEng','TotEng','KinEng']:
        thermo_property = thermo_property - min(thermo_property)
     
    start = bfp.get_nearestindex(ps, key='first')
    end = bfp.get_nearestindex(ps,key='last')
    
    plt.plot(ps[start:end],thermo_property[start:end],color='red')
    
    if showstep:
        plt.xlabel('Iteration')
    else:
        plt.xlabel(ne.getlab('time'))
        
    plt.ylabel(ne.getlab(plot_property))
    plt.title(directory[directory.find('Base_Oil'):]+'\n\n')
    # plt.text()
    # plt.xlim(0,1000+10)
    # plt.plot([700]*10,np.linspace(0,1.2,10),'--')
    # plt.savefig('python_outputs\\thermoplot\\thermoplot-'+ne.randstr(), dpi=300,bbox_inches='tight')
    plt.show()


