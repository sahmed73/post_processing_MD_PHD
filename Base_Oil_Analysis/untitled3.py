# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 13 17:00:10 2023
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def myfunction(arr):
    arr[0]=100
    
b = [1,2,3,4,5,6,7,8]
myfunction(b)
print(b)