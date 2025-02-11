# -*- coding: utf-8 -*-
"""
Created on Mon May 22 18:22:05 2023

@author: arup2
"""

import magnolia.formconverter as fc

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\A\Structure_A'

filename = r'\50A_300_O2_06_dense.xsd'

xsd_path = directory+filename

atom_type_order = ['H','C','O']
fc.xsd2data(xsd_path,atom_type_order)