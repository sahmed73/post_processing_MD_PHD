# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec 29 03:59:16 2023
"""

import magnolia.MD_Converter as mdc
import numpy as np

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250'

for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    filename  = 'bonds.reaxc'
    bondfile = directory+sim+'\\'+filename
    atomsymbols = ['H','C','O']
    
    for cutoff in np.array(range(35,80,5))/100:
        print(sim,cutoff)
        mdc.bond_to_speciecfile(bondfile, atomsymbols, cutoff)