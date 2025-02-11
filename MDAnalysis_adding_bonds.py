# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 19:11:36 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import MDAnalysis as mda

def mda_bond(bonded_dict, timestep):
    bonds = []
    tstep = timestep
    for atom in bonded_dict[tstep]:
        for neigh in bonded_dict[tstep][atom]:
            bonds.extend([(atom-1,neigh-1)])
    
    return bonds

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

neighbours,atomtypes = bfp.get_neighbours(bondfile)

#######################
trajectory_file = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\oxidation.lammpstrj'
topology_file = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\DataFile\50D_300_O2.data'
# Create a Universe object
u = mda.Universe(topology_file, trajectory_file,format='LAMMPSDUMP',dt=0.25)
traj = u.trajectory
#######################

for ts in traj:    # parsing through the dump file
    tstep = ts.data['step']    # saving the current timestep value
    
    u.add_TopologyAttr('bonds', mda_bond(neighbours, tstep))   # adding bond info to the topology
    ubonds = u.bonds    # Bond information of current step
#%%----
print(u.bonds)
for ub in traj.bonds:
    print(ub)