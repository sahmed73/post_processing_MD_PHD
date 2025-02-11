# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Sep  6 15:12:24 2024
"""
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import matplotlib.pyplot as plt
import numpy as np

topo = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\DataFile\Bulk\100_PAOr_50_A0001_density=0.2.data"
trej = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\REACTER\PAOr+Antioixdants\set001\A0001\Reaction\Sim-7_TSD=1.80\reaction.unwrapped.lammpstrj"

u = mda.Universe(topo, trej, format='LAMMPSDUMP')
MSD = msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
MSD.run()

msd =  MSD.results.timeseries
#%%


nframes = MSD.n_frames
timestep = 1  # This needs to be the actual time between frames
lagtimes = np.arange(nframes) * timestep/1000  # Make the lag-time axis


np.save("msd_A0001", np.column_stack((lagtimes, msd)))
fig, ax = plt.subplots(dpi=350)

# Plot the actual MSD
ax.plot(lagtimes, msd, color="black", linestyle="-", label=r'3D random walk')

exact = lagtimes * 6

# Plot the exact result
# ax.plot(lagtimes, exact, color="black", linestyle="--", label=r'$y=2 D\tau$')

ax.set_xlabel('Time (ns)')
plt.ylabel('MSD ($Å^2$)')
# Add legend and show the plot
ax.legend()
plt.show()
