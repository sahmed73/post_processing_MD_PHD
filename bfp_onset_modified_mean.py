# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 00:06:23 2023
@author: arup2
"""
import matplotlib.pyplot as plt
import magnolia.bondfile_parser as bfp
import magnolia.log_parser_FHB as lfp
import pandas as pd
import time
import numpy as np
from scipy.optimize import curve_fit
import random
import sys

start_time = time.time()
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\AO_Oxidation_with_Protocols\AO1\50_AO1_400_O2\Onset'

swsc_avg = pd.DataFrame()
timestep = 0.25
print('Timestep:',timestep)

for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
    d = directory+sim
    filename = '\\bonds.reaxc'
    bondfile_path   = d+filename
    bondpickle_path = d+'\\bonds.pickle'
    result = bfp.loadpickle_or_execute(bondpickle_path,bfp.get_neighbours,bondfile_path,bo=True,mtypes=True)
    #dropping bo and mtypes and adding atomsymbols
    atominfo = tuple(result[:2])+(['H','C','O'],) 

    neighbours,atomtypes,atomsymbols = atominfo
    
    species_pickle_path = d+'\\species.pickle'
    data = bfp.loadpickle_or_execute(species_pickle_path, bfp.get_SpeciesCountAtEveryTimestep,*atominfo,step2ps=True,step2psargs=[timestep])
    
    main_molecule = 'H12C10O3'#'H34C26O4'#
    swsc = pd.DataFrame(data).fillna(0).loc[main_molecule,:]
    print(swsc.head)
    swsc_avg = pd.concat([swsc_avg, swsc], axis=1)
#%%------------------------------------------------------
#swsc = pd.DataFrame(data).fillna(0).loc[,1250:2200]
swsc_avg = swsc_avg.fillna(0)
print(swsc_avg)
swsc = swsc_avg.mean(axis=1)
std = list(swsc_avg.std(axis=1).values)
print(swsc_avg.mean(axis=1))
#%%----------Extract Ramping Rate from logfile-----------
logfilename = '\\Sim-1\\log.lammps'
logfile     = directory+logfilename
thermo = lfp.thermo_dict(logfile, 1)
m,c = lfp.tempramp(thermo, timestep)

print('Ramping Rate (m):',round(m,2))
print('Initial Temperature (c):',round(c,2))
#%%----------------------------------
initial_molecule_count = 50

def fit_function(x, A, B, C, D):   
    y = C+((D-C)/(1+np.exp((x-A)/B)))
    return y

def inv_fit_function(y,A,B,C,D):
    x = A + B*np.log((D-y)/(y-C))
    return x



x = np.array(list(swsc.index))
y = np.array(list(swsc.values))
ramp = m
T0   = c
x = T0+ramp*x
print(len(x))
start = bfp.get_nearestindex(x, 0)
end = -1#bfp.get_nearestindex(x, x[-1])
x =x[start:end]
y = y[start:end]
std = std[start:end]
print(len(x))

################--Error Bars--###########################
# ex,ey,estd = [],[],[]
# for i in range(len(x)):
#     if i%100==0:
#         ex.append(x[i])
#         ey.append(y[i])
#         estd.append(std[i])
################--Individual Onset--#####################
onset_list=[]
for i in range(3):
    single = swsc_avg.iloc[:,i]
    xx = np.array(list(single.index))
    xx = T0+ramp*xx
    yy = np.array(list(single.values))
    yy[np.isnan(yy)]=0
    
    xx =xx[start:end]
    yy = yy[start:end]
    
    inital_guess = [2000,107,0,50]
    p, _ = curve_fit(fit_function, xx, yy,inital_guess,maxfev=2000)
    onset = inv_fit_function(49, *p)
    print('Onset-{}: {:.1f}'.format(i+1,onset))
    onset_list.append(onset)
onset_array = pd.DataFrame(onset_list)
print('\n-----Statistics-----\n',onset_array.describe(),'\n--------------------')
print(onset_array.describe().dtypes)
#################--Average Onset Curve Fit--###############
inital_guess = [2200,107,0,50]
parameters, covariance = curve_fit(fit_function, x, y,inital_guess,maxfev=2000)
print('A: {}\nB: {}\nC: {}\nD: {}\n'.format(*parameters))
fit_y = fit_function(x, *parameters)

factor = 0.98
whole  = 49#factor*initial_molecule_count
onset = inv_fit_function(whole, *parameters)
print('Onset: ',onset)


###########--Ploting--##############################
plt.scatter(x,y,marker='s', label='data',color='gray')
plt.plot(x, fit_y, '-', label='fit',color='r')
plt.plot([onset]*len(y),y,color='black')
#plt.errorbar(ex, ey, yerr=estd,fmt='o',color='black')
#plt.xticks(np.arange(1000, 2600, 500))
plt.xlabel('Temperature (K)')
plt.ylabel('Number of molecule')
plt.legend()


############-Plot Saving--############################
rand = random.randint(100000, 999999)

plt.title(directory[directory.find('AO_Oxidation'):]+'\nIndividual Onset: {:0.1f}, {:0.1f}, {:0.1f} \n{}\n A={:.0f},  B={:0.1f}, C={:0.1f}, D={:0.1f}\n\nOnset = {:0.1f} K'.format(*onset_list,onset_array.describe(),*parameters,onset))

sss = directory[directory.find('00ps')-1:].replace('\\', '-')+'-'
plt.savefig('python_outputs\\figures\\fit-onset-'+sss+str(rand),dpi=300,bbox_inches='tight')
plt.show()
print(np.array(bfp.step2picosecond(list(neighbours.keys()), timestep)))

runtime = time.time()-start_time
print('----Total run Time: {} min {:.0f} sec---------'.format(int(runtime/60),runtime%60))