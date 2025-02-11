# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Mar 14 22:03:26 2024
"""

import magnolia.MD_Converter as mdc
import pandas as pd



dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Material Structures'

pdbfile=dirr+'\\PAO_Radical.pdb'

atomtype_order = ['H','C','O']
box_size       = None # [0., 48., 0., 48., 0., 48.] # xlo, xhi, ylo, yhi, zlo, zhi
box_incr       = [0., 2., 0.,  2., 0.,  2.]

df = mdc.pdb2lmpdata(pdbfile, atomtype_order=atomtype_order,
                     box_size=box_size, box_incr=box_incr)