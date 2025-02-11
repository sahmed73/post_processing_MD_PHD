# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sun Sep 29 19:41:18 2024
"""

import magnolia.MD_Converter as mdc
import numpy as np
import pandas as pd
import magnolia.bondfile_parser as bfp
import magnolia.speciesfile_parser as sfp

dirr=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI-V2\1_ARLOI_24-09-12--00-35-24\20_PAO-OH_15_A0009\Equilibrate\Sim-1"
filename=r"\bonds.reaxc"
bondfile=dirr+filename

mdc.bond2speciecfile(bondfile, 'HCO', cutoff=0.30)