# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Feb 28 02:13:37 2024
"""
from magnolia.MD_Converter import *

###-User input-##################
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Material Structures'

filename = r'\20_PAO-OH_15_B.xsd'

xsdfile = directory+filename

atom_type_order = ['H','C','O']
xsd2lmpdata(xsdfile, atom_type_order)
# pass atom_type_order serially 
# If you do not pass the atom_type_order then the code will assign random type id for each element