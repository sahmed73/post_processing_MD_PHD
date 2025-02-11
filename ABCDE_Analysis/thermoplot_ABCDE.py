# -*- coding: utf-8 -*-
"""
Created on Mon May 22 23:33:13 2023

@author: arup2
"""

import magnolia.log_parser_FHB as lfp
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import numpy as np
import magnolia.needless_essential as ne
from scipy.optimize import curve_fit
plt.style.use('classic')

ABCDE = 'DBC'
fontsize = 15
extra = ['NPT_Only', 'NPT_Only','NPT_Only_06_dense']
marker = ['s', '^', 'v', 'o', '*', '>', '<']
for i, x in enumerate(ABCDE):
    
    directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\{}\{}_300_O2\Equilibrate\{}'.format(
        x, x, extra[i])
    filename = '\\log.lammps'
    logfile = directory+filename
    timestep = 0.25
    print('Molecule', x)
    thermo = lfp.thermo_dict(logfile, 2)

    plot_property = 'Density'

    steps = thermo['Step']
    ps = np.array(bfp.step2picosecond(steps, timestep))

    thermo_property = np.array(thermo[plot_property])
    print(plot_property, thermo_property[-1], '\n')
    if plot_property in ['PotEng', 'TotEng', 'KinEng']:
        thermo_property = thermo_property - min(thermo_property)

    start = bfp.get_nearestindex(ps, key='first')
    endpoint = 400
    end = bfp.get_nearestindex(ps, endpoint)
    
    if x == 'D': x='A'
    plt.plot(ps[start:end], thermo_property[start:end],
             label='Molecule '+x, marker=None, markevery=800)

plt.xlabel(ne.getlab('time'),fontsize=fontsize)
plt.xlim(-5, endpoint+5)
legend = plt.legend(loc="center right", edgecolor="black")
# legend.get_frame().set_alpha(None)
# legend.get_frame().set_facecolor((0, 0, 1, 0.1))
#plt.yaxis.set_ticks(np.arange(0, 85, step=20))
plt.ylabel(ne.getlab(plot_property),fontsize=fontsize)
plt.title('Density plot during NPT for Molecule A, B and C\n')
plt.savefig('..\\python_outputs\\thermoplot\\thermoplot-' +
            ne.randstr(), dpi=300, bbox_inches='tight')
