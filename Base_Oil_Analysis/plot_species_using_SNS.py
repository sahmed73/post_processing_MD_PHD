# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec 15 02:01:00 2023
"""

import seaborn as sns
import matplotlib.pyplot as plt
import magnolia.speciesfile_parser as sfp
import re
import pandas as pd

def get_element_count(formula, element):
    # Regular expression: element symbol followed by one or more digits
    pattern = f"{element}(\d+)"
    
    # Search for the pattern in the formula
    match = re.search(pattern, formula)
    
    # If the pattern is found, convert the number to an integer and return it
    if match:
        return int(match.group(1))
    else:
        # Return 0 if the element is not found or has no explicit number
        return 0



dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600"
specfilepath = dirr+'\\Sim-1\\species.out'

species = sfp.get_species_count(specfilepath).T

carbon = pd.DataFrame()
for i in range(100):
    req_col = []
    for col in species.columns:
        count = get_element_count(col,'C')
        if count==i:
            req_col.append(col)
    s = species[req_col].sum('columns')
    if s.sum()>0.0:
        carbon[f'C{i}']=s.copy()