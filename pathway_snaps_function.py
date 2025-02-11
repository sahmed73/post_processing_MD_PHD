# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 02:53:42 2023

@author: arup2
"""

import magnolia.bondfile_parser as bfp


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\B\B_300_O2\Production\1250\Sim-1'
filename  = '\\bonds.reaxc'
bondfilepath  = directory+filename

atomsymbols = ['H','C','O']
bonddata = bfp.parsebondfile(bondfilepath,ALL=True,pkl='yes')
#%%----------------------------------------------
neighbours = bonddata['neighbours']
bondorders = bonddata['bondorders']
atypes     = bonddata['atypes']
mtypes     = bonddata['mtypes']
charge     = bonddata['charge']
nlp        = bonddata['nlp']


seek = 'H2O'
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\Python_analysis\{}.txt'.format(seek)

result = bfp.species_to_molecule(bonddata,atomsymbols, seek, source=True, dump = savedir,shortinfo=True)