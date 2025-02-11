# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Feb 20 09:14:10 2024
"""
import numpy as np

u  = [1,2,3]
v  = [-1,0,4]

# distance
r  = np.sqrt(np.sum((u-v)**2))
r2  = np.sqrt((u[0]-v[0])**2 + (u[1]-v[1])**2 + (u[2]-v[2])**2)

print(r,r2)