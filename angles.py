# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 03:43:16 2023

@author: arup2
"""
import magnolia.dumpfile_parser as dfp
import magnolia.bondfile_parser as bfp

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\D\D_300_O2\Production\1200K'
dumpfilepath = directory+'\\'+'oxidation.lammpstrj'
bondfilepath = directory+'\\'+'bonds.reaxc'

dumpdata = dfp.parsedumpfile(dumpfilepath)
# bonddata = bfp.parsebondfile(bondfilepath)

atoms = {2568, 2566, 2567}
positions = dumpdata['position']
atomtypes = dumpdata['atypes']

atoms = [atomtypes[x] for x in atoms]
print(atoms)

for step, pos in positions.items():
    pass