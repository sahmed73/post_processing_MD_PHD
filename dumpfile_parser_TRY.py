# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Mar 31 05:49:56 2024
"""
import magnolia.dumpfile_parser as dfp
import pandas as pd

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO_Radical\20_PAO_Radical_20_BHT_Soria\Equilibrate\100NVT+400NPT'

filename = '\\equilibrated.nvt.lammpstrj'
dumpfile = dirr+filename

dumpdata = dfp.parsedumpfile(dumpfile)
