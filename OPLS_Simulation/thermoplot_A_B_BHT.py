# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Apr  7 14:57:07 2024
"""
import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
mplt.custom_plot_features(minorticks=True)


antioxidants = ['A','B','BHT']
for ao in antioxidants:
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\OPLSAA\PAO_Antioxidant\PAO+{}\10PAO+3{}\Equilibrate\Energy'.format(ao,ao)
    filename = '\\log.lammps'
    logfile = dirr+filename
    
    thermo = lfp.thermo_panda(logfile, serial=2, timestep=0.25)
    print(thermo.head())
    
    prop_1 = 'Time'
    prop_2 = 'PotEng'
    x, y = thermo[prop_1], thermo[prop_2]
    
    prop_list = ['PotEng','TotEng','KinEng', 'Time']
    if prop_1 in prop_list:  x = x  - x.min()
    if prop_2 in prop_list:  y = y  - y.min()
    plt.plot(x,y, label=ao)
    # plt.xlim(right=1000)
plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel(f'{ne.getlab(prop_2)}')
plt.legend(loc='center right')
mplt.saveplot(folder='thermoplot',name='thermo')