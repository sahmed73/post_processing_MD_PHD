# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Apr  3 10:25:39 2024
"""

import numpy as np

n = np.loadtxt('hudai.txt',dtype=int)

for i in range(1,92+1):
    if i in n: continue
    print(i)