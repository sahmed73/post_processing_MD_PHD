# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:50:32 2023

@author: arup2
"""
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import pandas as pd

#x =      [20,    10,     4,       2,       1,     0.5][0:]
#y=       [1990.2,1964.9, 1779.2,  1909.1,  1699.3,1490][0:]
#y =      [2164,  1948.1, 1807,    1837.9,  1736.1,  1490]#[0:]
error =  [117,   98.4,   85,      67.4,    51.3,  52.2][0:]
#error =  [74.5,  56.7, 54.5,    35.5,  49.4,  52.2]#[0:-1]

#onset_sameX = {20:1990, 10:1970.5,4:1779, 2:1909, 1:1699, 0.5:1490}
onset_sameX_1000_2500 = {1:1754.8,2:1800,4:1789,10:1964,20:2025}
error_sameX_1000_2500 = {1:49,2:28.6,4:42,10:4.3,20:50}
x = list(onset_sameX_1000_2500.keys())
y = list(onset_sameX_1000_2500.values())
error = list(error_sameX_1000_2500.values())

for a,b in zip(x,y):
    print(a,'=',b)
#computational_time = ['10h 6m 58s/100','62h 26m 15s/459','29h 42m 53s/750','72h 0m 0s/1812']
#computational_time = [100/10.1,459/62.4,750/29.7,1812/72]
#computational_time = np.array(list(onset_sameX.keys()))
#y = 2/computational_time


ran = random.randint(10000,99999)
savefile = 'ramp_vs_onset'+'-'+str(ran)
plt.semilogx(x,y,color='black',marker='s')
plt.errorbar(x, y, yerr=error,fmt='.',color='black')
#plt.yticks(np.arange(0, 5, step=1))
#plt.title('Ramping temperature from 1000K to 3000K',fontsize=8)
plt.xlabel('Ramping rate (K/ps)')
plt.ylabel('Onset (K)')
plt.ylabel('Computational time (ns)')
plt.savefig('figures\\'+savefile,dpi=300,bbox_inches='tight')