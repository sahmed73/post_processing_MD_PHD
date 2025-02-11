# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Jan 25 04:10:02 2024
"""

import magnolia.MD_Converter as mdc
import numpy as np

base_oil = 'Squalane'
directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600".format(base_oil,base_oil)

for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    filename  = 'bonds.reaxc'
    bondfile = directory+sim+'\\'+filename
    atomsymbols = ['H','C','O']
    
    for cutoff in np.array(range(35,80,5))/100:
        print(sim,cutoff)
        mdc.bond_to_speciecfile(bondfile, atomsymbols, cutoff)