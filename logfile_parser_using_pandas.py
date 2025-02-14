# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  6 15:26:03 2024
"""
import magnolia.log_parser_SA as lfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne
import magnolia.plot_template as mplt
import numpy as np

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=18


dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0004\Eq\Sim-1_Interaction_Energy+MSD"
filename = r'\log.lammps'

logfile = dirr+filename

thermo = lfp.thermo_panda(logfile, serial=3)

prop_1 = 'Time'
prop_2 = 'Density'#'c_PAO_AO_interaction'
x, y = thermo[prop_1], thermo[prop_2]

start = 0
end   = None

x, y = x.iloc[start:end], y.iloc[start:end]

fig, ax = plt.subplots(dpi=300)
ax.plot(x,y,c='tab:red')
plt.xlabel(f'{ne.getlab(prop_1)}')
plt.ylabel(f'{ne.getlab(prop_2)}')
# plt.xlabel('Time (ns)')
# plt.ylabel('Number of reactions')