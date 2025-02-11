# -*- coding: utf-8 -*-
"""
Created on Mon May 15 23:07:14 2023

@author: arup2
"""

import MDAnalysis as mda
import sys
import math
import matplotlib.pyplot as plt

# Load the trajectory and topology files

def compute_angle(p,q,r):
    A  = [x-y for x,y in zip(p,q)]
    B = [x-y for x,y in zip(p,r)]
    
    AB = sum([x*y for x,y in zip(A,B)])
    mA = sum([x**2 for x in A])**0.5
    mB = sum([x**2 for x in B])**0.5
    
    cos = AB/(mA*mB)
    angle = 180*math.acos(cos)/(2*math.pi)
    return angle

trajectory_file = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K\oxidation.lammpstrj'

topology_file = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\DataFile\50D_300_O2.data'

# Create a Universe object
u = mda.Universe(topology_file, trajectory_file,format='LAMMPSDUMP',dt=0.25)

atoms = [2822]#, 1526, 1527] # O=C=O
criteria = "id "+" ".join(map(str,atoms))
# print(criteria)
selected_atoms = u.select_atoms(criteria)
# print(selected_atoms)

angles = []
for ts in u.trajectory:
    positions = selected_atoms.positions
    # print(selected_atoms.types)
    temp = []
    for atom in selected_atoms:
        print(ts.data['step'])
        print(atom.position)
        if atom.type == '2':
            p = atom.position
        else:
            temp.append(atom.position)
    # q, r = temp
    # angle = compute_angle(p,q,r)
    # angles.append(angle)
    # print('{:0.2f}'.format(angle))

# print(sum(angles)/len(angles))    
# xx = range(len(angles))
# plt.plot(xx,angles)
    