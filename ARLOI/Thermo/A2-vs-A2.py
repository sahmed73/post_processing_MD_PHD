# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Oct 22 01:46:52 2024
"""


import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import pandas as pd

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=15

fig, ax = plt.subplots(dpi=350)

#---NPT1:600K,2atm(100ps)
dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\7_A0001-5_eq_NPT-NPT-NVT_Pro_2ramp\A0002\Equilibration"
logfile = dirr+r'\log.lammps'
thermo1 = lfp.thermo_panda(logfile, serial=2, timestep=0.25, zero_ref='time')
#---NPT2:600->300K,2->1atm(50ps)---NPT3:300K,1atm(100ps)
dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\8_ARLOI_Down_to_300K\A0002\Equilibration\Sim-1"
logfile = dirr+r'\log.lammps'
thermo2 = lfp.thermo_panda(logfile, serial=[1,2], timestep=0.25)

thermo=pd.concat([thermo1, thermo2], axis='rows')
ax.plot(thermo['Time'],thermo['Density'], label='Eq$_{new}$')

#---NPT:300K,1atm(400ps)
dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\4_ARLOI_24-09-07--22-34-03\20_PAO-OH_15_A0001\Equilibrate\Sim-1"
logfile = dirr+r'\log.lammps'
thermo = lfp.thermo_panda(logfile, serial=3, timestep=0.25, zero_ref='time')
ax.plot(thermo['Time'],thermo['Density'], label='Eq$_{old}$')

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Density (g/cm$^3$)')
plt.legend()
plt.grid()