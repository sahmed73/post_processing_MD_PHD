# -*- coding: utf-8 -*-
"""
Created on Wed May 31 23:28:53 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys

bf1 = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\100_Less_TempRamp\bonds.reaxc'
bf2 = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\100_Less_TempRamp_Continuation\bonds.reaxc'
mf  = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\merge_bonds.reaxc'

# path = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B'
# bf1 = path+r'\bonds_1.reaxc'
# bf2 = path+r'\bonds_2.reaxc'
# bf3 = path+r'\bonds_3.reaxc'
# mf  = path+r'\marge_bonds.reaxc'
# real = path+r'\bonds.reaxc'

bfp.merge_bondfiles(bf1,0,4498000,bf2,4499000,'last',file=mf)
neighbours, at = bfp.get_neighbours(mf)
#%%
steps = list(neighbours.keys())
ps = bfp.step2picosecond(steps, 0.25 )
index = range(len(steps))
plt.plot(index,ps)
# r = bfp.text_isequal(real, mf,strip=True)
# print(r)
