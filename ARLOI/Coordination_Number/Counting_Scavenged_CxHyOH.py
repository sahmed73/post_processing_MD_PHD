# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Oct 23 16:10:42 2024
"""

'''
A function to identify and count CxHyOH molecules,
marking any carbon chain with an OH group as a scavenged radical.
'''

import magnolia.bondfile_parser as bfp
import pandas as pd
import matplotlib.pyplot as plt

### Directory ###
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Production\fixed-300K"
filename  = "\\bonds.out"
bondfile = directory+filename
atominfo = bfp.parsebondfile(bondfile, mtypes=True)
#%%
def count_scavenged(atominfo):
    neighbours = atominfo['neighbours']
    atypes  = atominfo['atypes']
    mtypes = atominfo['mtypes']
    
    # types
    hydrogen_type=1
    carbon_type=2
    oxygen_type=3
    
    scavenged = {}
    
    scavenged_with_AO_hydroxyl_h={}
    scavenged_with_other_h={}
    
    
    # is the H from AO hydroxyl
    is_from_AO_hydroxyl={} # bool
    for step, neigh in neighbours.items():
        for parent, children in neigh.items():
            is_from_AO_hydroxyl[parent]=False
            
            if atypes[parent]==hydrogen_type and mtypes[parent]>init_num_of_radicals:
                child = children[0]
                if atypes[child]==oxygen_type:
                    is_from_AO_hydroxyl[parent]=True
                if len(children)>1:
                    print("Warning: some H have more than one child")
        break # only checking in the first step
                
    for step, neigh in neighbours.items():
        scavenged[step]=0
        scavenged_with_AO_hydroxyl_h[step]=0
        scavenged_with_other_h[step]=0
        
        for parent, children in neigh.items():
            # count if it lookes like -C-O-H
            if atypes[parent]==oxygen_type and mtypes[parent]<=init_num_of_radicals:
                children_types=sorted([atypes[x] for x in children])
                if children_types==[hydrogen_type,carbon_type]: # H and C
                    scavenged[step]+=1
                    
                    for child in children:
                        if atypes[child]==hydrogen_type:
                            if is_from_AO_hydroxyl[child] is True:
                                scavenged_with_AO_hydroxyl_h[step]+=1
                            else:
                                scavenged_with_other_h[step]+=1
    
    return pd.Series(scavenged), pd.Series(scavenged_with_AO_hydroxyl_h), pd.Series(scavenged_with_other_h)


timestep  = 0.25
ramp_rate = 4
initial_temp = 300
init_num_of_radicals = 20
init_num_of_AOs = 15
 
scavenged, hydroxyl, others = count_scavenged(atominfo)
time=scavenged.index*timestep/1000

scavenged.index=time
hydroxyl.index=time
others.index=time

plt.rcParams['font.size']=16
fig, ax = plt.subplots(dpi=350)
ax.plot(scavenged, label='Total')
ax.plot(hydroxyl, label='Scavenged by AO-Hydroxyl H')
ax.plot(others, label='Scavenged by Other H')
ax.set_xlabel("Time (ps)")
ax.set_ylabel("Number of scavenged radicals")
ax.legend(fontsize=11,loc='lower right')

# Add secondary x-axis for temperature
def temp_function(x):
    """Convert time (ps) to temperature."""
    return initial_temp + x * ramp_rate

def time_function(x):
    """Convert temperature back to time (ps)."""
    return (x - initial_temp) / ramp_rate

secax = ax.secondary_xaxis('top', functions=(temp_function, time_function))
secax.set_xlabel('Temperature (K)')
#%%