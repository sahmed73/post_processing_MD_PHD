# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Jun 19 14:48:59 2024
"""

import magnolia.speciesfile_parser as sfp
import matplotlib.pyplot as plt

dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO-OH\20_PAO-OH_15_A_v2\Onset\Sim-1'
speciesfile = dirr+'\\species.out'

species = sfp.get_species_count(speciesfile, timestep=0.25)