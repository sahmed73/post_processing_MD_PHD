# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 12 23:13:41 2023
"""

import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import sys
#%%
baseoils = ['PAO4']
location = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_5A_Soria\Onset'
directory_list = [location.format(oil) for oil in baseoils]


## Global Variables
timestep        = 0.25
ramp_rate       = 4
initial_temp    = 300
baseoil_formula = {'PAO4': 'H62C30', 'Squalane': 'H62C30'}
imc             = 25

## Printing some variables
print('Ramp Rate',ramp_rate,'K/ps')
print('Initial Temperature:',initial_temp)

light_colors = ['#FF7779','#89CFF0','#CEFAD0']
dark_colors  = ['']


## Variable Initialization
mean_temp_list = []
mean_sc_list = []
mean_onset_list = []
sc_fit_list = []
error_list = []


for index, directory in enumerate(directory_list):
    baseoil = baseoils[index]
    
    onsets    = []
    tempss    = []
    scss      = []

    for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
        print(baseoil)
        bondfilepath = directory+sim+'\\bonds.reaxc'
        if not os.path.exists(bondfilepath):
            sys.exit('No file exists!!')
            
        bonddata = bfp.parsebondfile(bondfilepath)
        
        neighbours = bonddata['neighbours']
        atypes     = bonddata['atypes']
        asyms      = ['H','C','O']
        
        main_molecule = baseoil_formula[baseoil]
        
        temp,sc = bfp.get_speciesVStemp(main_molecule, neighbours, atypes,
                                        asyms,it=initial_temp,
                                        ramp=ramp_rate,ts=timestep)
        
        onset,_   = bfp.get_onset(temp, sc,imc=imc,ig=[1700,107,0,50])
        tempss.append(temp)
        scss.append(sc)
        onsets.append(onset)
        print('Onset:',onset,'\n')
  
    #### Addition:Plotting Shaded Region for 3 sets of data ####
    max_data = []
    min_data = []
    for j in range(len(scss[0])):
        max_data.append(max(scss[0][j],scss[1][j],scss[2][j]))
        min_data.append(min(scss[0][j],scss[1][j],scss[2][j]))
    plt.fill_between(tempss[0], min_data, max_data, alpha=0.3,
                     color=light_colors[index])
    #########################################

    mean_temp         = sum(tempss)/len(tempss)
    mean_sc           = sum(scss)/len(scss)
    mean_onset,sc_fit = bfp.get_onset(mean_temp, mean_sc, imc=imc,ig=[1700,107,0,25])
    print('Onsets: {:.0f}, {:.0f}, {:.0f}'.format(*onsets))
    print('Mean Onset: {:.0f}'.format(mean_onset))
    
    #finding error
    error = 0
    for onset in onsets:
        error += (onset-mean_onset)**2
    error = (error/len(onsets))**0.5
    print('error: {:.1f}'.format(error))
    
    mean_temp_list.append(mean_temp)
    mean_sc_list.append(mean_sc)
    mean_onset_list.append(mean_onset)
    sc_fit_list.append(sc_fit)
    error_list.append(error)
    
#%%--------plotting-----
plt.style.use('classic')
color = ['r','b','g']
light_color = ['#FFDBE9','#CAE9F5','#CEFAD0']
label = ['Molecule B', 'Molecule A', 'Molecule B']
dashedTop = [imc-0.4]*3

n = len(directory_list)

for i in range(1+n):
    if i==0: continue
    i = i%n
    # adjust shift ##
    mean_sc_list[i] = mean_sc_list[i]*(imc/mean_sc_list[i][0])
    sc_fit_list[i] = sc_fit_list[i]*(imc/sc_fit_list[i][0])
    #################
    
    # plt.scatter(mean_temp_list[i],mean_sc_list[i],
    #             marker='s',color=light_color[i],s=50,alpha=1)
    plt.plot(mean_temp_list[i],sc_fit_list[i], 
             label=label[i],color=color[i],linewidth=2)
    plt.plot([mean_onset_list[i]]*1000,np.linspace(-5,dashedTop[i],1000),'--',
             color=color[i],linewidth=2)
    
plt.ylim(-1,26)
plt.xlim(00,3550)
plt.xlabel('Temperature (K)',fontsize=15)
plt.ylabel('Number of molecules',fontsize=15)
plt.legend()
#plt.grid(visible=True)
plt.savefig('..\\python_outputs\\onset\\onset-'+str(random.randint(1000000,9999999)), dpi=1000,bbox_inches='tight')
#-------------------