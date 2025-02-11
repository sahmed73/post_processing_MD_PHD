# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Mar 31 08:48:51 2024
"""

import pandas as pd


dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO_Radical\Eng_min_forced_dens_20_PAO_Radical_20_BHT_Soria\DataFile'
filename = '\\eng_min_forced_dens_20_PAO_Radical_20_BHT.data'

datafile = dirr+filename

data = pd.read_csv(datafile,delimiter=r'\s+',
                   skiprows=17, header=None).iloc[:,:-3]

adata = data.copy() - [0,0,0,0,12,0,0]

file_path = dirr+'output.data'
with open(file_path, 'w') as file:
    file.write(adata.to_string())