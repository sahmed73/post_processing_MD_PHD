# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 04:33:13 2023

@author: shihab
"""

import matplotlib.pyplot as plt
import magnolia.bondfile_parser as bfp
import magnolia.log_parser_FHB as lfp
import pandas as pd
import time
import sys
import numpy as np
from scipy.optimize import curve_fit
import random

start_time = time.time()
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_200_O2_Cubic_68\Onset\Sim-1'

filename = '\\bonds.reaxc'
bondfile  = directory+filename
neighbours, atomtypes = bfp.get_neighbours(bondfile,atypes=True)

data = bfp.get_SpeciesCountAtEveryTimestep(neighbours, atomtypes, 'HCO',step2ps=True,step2psargs=[0.25])
#%%----------Extract Ramping Rate from logfile-----------
logfilename = '\\log.lammps'
logfile     = directory+logfilename
thermo = lfp.thermo_dict(logfile, 1)

def ramp_fit(x,m,c):
    y = m*x+c
    return y

steps = thermo['Step']
ps = np.array(bfp.step2picosecond(steps, 0.25))
temp = thermo['Temp']
x = np.array(ps)
y = np.array(temp)
p, cov = curve_fit(ramp_fit, x, y)

m = p[0]
c = p[1]
print('Ramping Rate (m):',round(m,2))
print('Initial Temperature (c):',round(c,2))

#%%-----------------------------------

def fit_function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_fit_function(y,A,B,C,D):
    x = A + B*np.log((D-y)/(y-C))
    return x

def second_derivative(x,A,B):
    imc = 50 #initial molecule count
    ex = np.exp((x-A)/B)
    ddy1 = 2*ex**2/((1+ex)**3)
    ddy2 = ex/((1+ex)**2)
    
    ddy = imc*(ddy1-ddy2)
    return ddy
#%%------------------------------------------------------
main_molecule = 'H12C10O3'
swsc = pd.DataFrame(data).fillna(0).loc[main_molecule,:]

x = np.array(list(swsc.index))
y = np.array(list(swsc.values))
ramp = m
T0   = c
x = T0+ramp*x

inital_guess = [2200,107,0,50]
parameters, covariance = curve_fit(fit_function, x, y,inital_guess,maxfev=2000)
fit_A = parameters[0]
fit_B = parameters[1]
fit_C = parameters[2]
fit_D = parameters[3]
print('A: ',fit_A)
print('B: ',fit_B)
print('C: ',fit_C)
print('D: ',fit_D)
fit_y = fit_function(x, *parameters)


#fit_2y= second_derivative(x, fit_A, fit_B)
#print("Max 2nd derivative: ",min(fit_2y))
#print("Max Corresponding x: ",inv_fit_function(-min(fit_2y), fit_A, fit_B))

#####################################
find_something = fit_A
closest_id = min(range(len(x)), key=lambda i: abs(x[i]-find_something))
tol = 0.98*50 #2% tolarance

onset = inv_fit_function(49, *parameters)
#####################################

plt.scatter(x,y,marker='s', label='data',color='gray')
plt.plot(x, fit_y, '-',color='r', label='fit')
plt.plot([onset]*len(y),y,color='Black')
plt.xticks([v for v in range(1000,3100+1,500)])
plt.xlabel('Temperature (K)')
plt.ylabel('Number of molecule')
plt.legend()

func = '$y=\\frac{{50}}{{1+exp(\\frac{{x-A}}{B})}}^m$'
rand = random.randint(100000, 999999)

plt.title(directory[directory.find('AO_Oxidation'):]+'\nFit function: {}\n A={:.0f},  B={:0.1f}, C={:0.1f}, D={:0.1f}\n\nOnset = {:0.1f} K'.format(func,*parameters,onset))

sss = directory[directory.find('00ps')-1:].replace('\\', '-')+'-'
plt.savefig('python_outputs\\figures\\fit-onset-'+sss+str(rand),dpi=300,bbox_inches='tight')
#plt.grid()
plt.show()



print('Onset: ',onset)


runtime = time.time()-start_time
print('----Total run Time: {} min {:.0f} sec---------'.format(int(runtime/60),runtime%60))
