# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Mar 26 03:12:48 2024
"""

import pandas as pd

dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing"

df = pd.read_csv(dirr+'\\charge_csv.csv',header=None)
charges = df.iloc[:,3]
