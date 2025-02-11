# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 14 12:34:17 2024
"""
import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt

# mplt.custom_plot_features(minorticks=True)
plt.style.use('default')
plt.rcParams['font.size']=15

fig, ax = plt.subplots(dpi=350)

#---NPT1:600K,2atm(100ps)---NPT2:600->300K,2->1atm(50ps)---NPT3:300K,1atm(100ps)
dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\8_ARLOI_Down_to_300K\A0001\Equilibration\Sim-2_48-cores"
logfile = dirr+r'\log.lammps'
thermo = lfp.thermo_panda(logfile, serial=[2,3,4], timestep=0.25, zero_ref='time')
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