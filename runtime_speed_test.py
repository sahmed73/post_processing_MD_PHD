# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:51:23 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp
import magnolia.needless_essential as ne
import time


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_200_O2\Production\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

start_time = time.time()
neighbours_v1,atomtypes_v1 = bfp.get_neighbours(bondfile,atypes=True)
runtime_v1 = time.time()-start_time
print('Version 1:')
ne.print_runtime(runtime_v1)

start_time = time.time()
neighbours_v2,atomtypes_v2 = bfp.get_neighbours_v2(bondfile)
runtime_v2 = time.time()-start_time
print('Version 2:')
ne.print_runtime(runtime_v2)
#%%--
neighbours_v2,atomtypes_v2 = bfp.get_neighbours(bondfile)
#print('speed gain by version 2:')
#print((runtime_v1)/runtime_v2,'times')