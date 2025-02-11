# -*- coding: utf-8 -*-
"""
Created on Mon May 22 19:03:16 2023

@author: arup2
"""

import magnolia.formconverter as fc
import numpy as np


for x in 'C':
    print('-'*20,x,'-'*20)
    
    data_path=r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\{}\Structure_{}\50{}_300_O2.data'.format(x,x,x)
    
    dd = fc.distance_between_all_pair(data_path)
    d = np.array(dd)
    #d = sorted(d)
    
    print('Min:',d.min())
    print('Max:',d.max())
    print('Mean:',d.mean())
    print('Std:',d.std())