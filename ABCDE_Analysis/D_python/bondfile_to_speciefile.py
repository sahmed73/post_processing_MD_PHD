# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec 29 03:59:16 2023
"""

import magnolia.MD_Converter as mdc
import numpy as np

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\Pure_50B\Production\1050'

for sim in ['\\Sim-1']: #,'\\Sim-2','\\Sim-3']:
    filename  = 'bonds.reaxc'
    bondfile = directory+sim+'\\'+filename
    atomsymbols = ['H','C','O']
    
    for cutoff in np.array(range(35,55,5))/100:
        print(sim,cutoff)
        mdc.bond2speciecfile(bondfile, atomsymbols, cutoff)