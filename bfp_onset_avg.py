# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 00:06:23 2023

@author: arup2
"""
import matplotlib.pyplot as plt
import magnolia.bondfile_parser as bfp
import pandas as pd
import time
import sys
import winsound

start_time = time.time()
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation\AO-1\50_AO1_100_O2\CHO_ReaxFF-2016\Onset\20000K-per-ns\Sim-'

swsc_avg = []
for case in [1,2,3]:
    d = directory+str(case)
    filename = '\\bonds.reaxc'
    bondfile  = d+filename
    neighbours, atomtypes = bfp.get_neighbours(bondfile,atypes=True)
    
    data = bfp.get_SpeciesCountAtEveryTimestep(neighbours, atomtypes, 'HCO',step2ps=True,step2psargs=[0.1])
    
    main_molecule = 'H12C10O3'
    swsc = pd.DataFrame(data).fillna(0).loc[main_molecule,:]
    print(swsc.head)
    swsc_avg.append(swsc)
#%%-------------------------------------------
swsc = ((swsc_avg[0]+swsc_avg[1]+swsc_avg[2])/3).fillna(0)
print(swsc)
print(swsc_avg[0].iloc[60],swsc_avg[1].iloc[60],swsc_avg[2].iloc[60],swsc.iloc[60])
#%%-----------------------------------
import numpy as np
from scipy.optimize import curve_fit
import random

def Gauss(x, A, B):
    y = 50/(1+np.exp((x-A)/B))
    return y

x = list(swsc.index)
y = list(swsc.values)

#####################Taking Few data################
xx = []
yy = []
for x_value,y_value in zip(x,y):
    x_value,y_value = round(x_value,1),round(y_value,1)
    t = 100 ###########################
    ramp = 20
    x_value = 1000+ramp*x_value
    xx.append(x_value)
    yy.append(y_value)
x = np.array(xx)
y = np.array(yy)
####################################################

parameters, covariance = curve_fit(Gauss, x, y,[2129,107])
fit_A = parameters[0]
fit_B = parameters[1]
print('A: ',fit_A)
print('B: ',fit_B)
fit_y = Gauss(x, fit_A, fit_B)
#####################################
find_something = fit_A
closest_id = min(range(len(x)), key=lambda i: abs(x[i]-find_something))
tol = 0.02*50 #2% tolarance
for i in range(1000000):
    right_id = closest_id + i
    left_id = closest_id - i
    if abs(y[right_id])<=tol  and (y[left_id]-50)<=tol:
        break
#####################################
onset = 2164#x[left_id]

plt.scatter(x,y,marker='s', label='data',color='gray')
plt.plot(x, fit_y, '-', label='fit',color='r')
plt.plot([onset]*len(y),y,color='black')
#plt.plot([x[right_id]]*len(y),y,color='black')
#plt.plot([x[closest_id]]*len(y),y,color='Blue')
plt.xticks([v for v in range(1000,3100+1,500)])
plt.xlabel('Temperature (K)')
plt.ylabel('Number of molecule')
plt.legend()

print('Onset: ',onset)
func = '$y=\\frac{{50}}{{1+exp(\\frac{{x-A}}{B})}}^m$'#'$y=e^{{Ax^{{B}}}}$'
rand = random.randint(100000, 999999)

plt.title(directory[directory.find('AO_Oxidation'):]+'\nFit function: {}\n A={:.0f},  B={:0.1f}\n\nOnset = {} K'.format(func,fit_A,fit_B,onset))

sss = directory[directory.find('00ps')-1:].replace('\\', '-')+'-'
plt.savefig('figures\\fit-onset-'+sss+str(rand),dpi=300,bbox_inches='tight')
#plt.grid()
plt.show()


runtime = time.time()-start_time
print('----Total run Time: {} min {:.0f} sec---------'.format(int(runtime/60),runtime%60))
sys.exit()
################### Beep after finish running #########################
winsound.Beep(1200, 100)
winsound.Beep(2000, 100)
winsound.Beep(1200, 100)