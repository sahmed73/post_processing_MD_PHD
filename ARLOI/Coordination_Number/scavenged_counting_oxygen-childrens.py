# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 21 14:57:19 2024
"""
import magnolia.bondfile_parser as bfp
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

### Directory ###
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\2_ARLOI_24-08-29--00-27-52\20_PAO-OH_15_A0002\Production\Sim-1"
filename  = "\\bonds.reaxc"
bondfile = directory+filename

atominfo = bfp.parsebondfile(bondfile, mtypes=True)
#%% count the coordination numbers
neighbours = atominfo['neighbours']
atypes  = atominfo['atypes']
mtypes = atominfo['mtypes']

timestep  = 0.25
ramp_rate = 4
initial_temp = 300

scavenged={}
from_AO_hydroxyl={}
from_others={}

sk_len=set()
for step, neigh in neighbours.items():
    time=step*timestep/1000
    temp=time*ramp_rate+initial_temp
    scavenged[temp]=0
    
    from_AO_hydroxyl[temp]=0
    from_others[temp]=0
    
    molecules=bfp.get_molecules(neigh)
    for parent, children in neigh.items():
        
        if mtypes[parent]<=20 and atypes[parent]==3:
            
            for molecule in molecules:
                if parent in molecule:
                    skeleton_len=[atypes[x] for x in molecule].count(2)
                    break
            
            
            coord_num=len(children)
            if coord_num==2:
                for child in children:
                    if atypes[child]==1 and mtypes[child]>20:
                        scavenged[temp]+=1
                        # print(skeleton_len)
                        sk_len.add(skeleton_len)
                        from_AO_hydroxyl[temp]+=1
                    elif atypes[child]==1:
                        # print(skeleton_len)
                        sk_len.add(skeleton_len)
                        scavenged[temp]+=1
                        from_others[temp]+=1

df=pd.Series(scavenged)
print(sk_len)
#%% plotting
fig, ax = plt.subplots(dpi=350)
ax.plot(df, label=f'$TRR={ramp_rate} K\cdot ps^{-1}$')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Counts')
ax.legend()
#%% count the H comes from antioxidants
df1=pd.Series(from_AO_hydroxyl)/2
df2=pd.Series(from_others)/2
fig, ax = plt.subplots(dpi=350)
ax.plot(df1,label='from_AO_hydroxyl')
ax.plot(df2,label='from_others')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('H counts')
ax.legend()