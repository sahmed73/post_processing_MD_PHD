# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 02:33:53 2023

@author: shihab
"""
import magnolia.bondfile_parser as bfp
import numpy as np
import matplotlib.pyplot as plt
import os
import random
#%%
lignins = {'A':'H12C10O3',
           'B':'H34C26O4',
           'C':'H44C29O2',
           'D':'H30C19O3',
           'E':'H14C11O4'}

pre = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE'
directory_list = [pre+r'\D\D_300_O2\Onset']#,
#                  pre+r'\B\B_300_O2\Onset',]
                  #pre+r'\C\C_300_O2\Onset_v2']
light_color = ['#FF7779','#89CFF0']
mean_temp_list = []
mean_sc_list = []
mean_onset_list = []
sc_fit_list = []
error_list = []
upper = []
lower = []

for i in range(len(directory_list)):
    directory = directory_list[i]

    timestep = 0.25
    index = directory.find('ABCDE')
    x = directory[index+6:index+7]
    ramp_rate      = 4
    initial_temp   = 300
    print('Molecule:',x)
    print('Ramp Rate',ramp_rate,'K/ps')
    print('Initial Temperature:',initial_temp)

    onsets    = []
    tempss    = []
    scss      = []

    for sim in ['\\Sim-1','\\Sim-2','\\Sim-3']:
        print('------------',sim[1:],'------------')
        d = directory+sim
        filename = '\\bonds.reaxc'
        bfpath   = d+filename # bond file path
        if not os.path.exists(bfpath):
            continue
        atominfo = bfp.get_neighbours(bfpath)
        neighbours,atomtypes = atominfo
        atomsymbols = ['H','C','O']
        
        main_molecule = lignins[x]
        print(main_molecule)
        
        temp,sc = bfp.get_speciesVStemp(main_molecule, *atominfo,
                                        atomsymbols,it=initial_temp,
                                        ramp=ramp_rate,ts=timestep)
        onset,_   = bfp.get_onset(temp, sc,imc=50,ig=[1700,107,0,50])
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
    upper.append(max_data)
    lower.append(min_data)
    #########################################

    mean_temp         = sum(tempss)/len(tempss)
    mean_sc           = sum(scss)/len(scss)
    mean_onset,sc_fit = bfp.get_onset(mean_temp, mean_sc, imc=50,ig=[1700,107,0,50])
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
import sys
plt.style.use('default')
plt.rcParams['font.size']=20
color = ['r','b','g']
label = ['Molecule A', 'Molecule B']
dashedTop = [49.55,49.55,49.55]
label_fontsize = 16

n = len(directory_list)

fig, ax = plt.subplots(dpi=300)
for i in range(0,n):
    # adjust shift ##
    mean_sc_list[i] = mean_sc_list[i]*(50/mean_sc_list[i][0])
    sc_fit_list[i] = sc_fit_list[i]*(50/sc_fit_list[i][0])
    #################
    
    plt.plot(mean_temp_list[i],sc_fit_list[i], 
             label=label[i],c=color[i],linewidth=2)
    plt.plot([mean_onset_list[i]]*1000,np.linspace(-5,dashedTop[i],1000),'--',
             c=color[i],linewidth=2)
    plt.fill_between(tempss[0], lower[i], upper[i], alpha=0.3,
                     color=light_color[i])

# plt.text(1220,15,'1351 K',color='b',rotation=90,fontsize=label_fontsize)
ax.text(925,15,'1104 K', color='r', rotation=90)
plt.ylim(-1,53)
plt.xlim(0,3550)
plt.xticks(range(0, 3550, 1000))
plt.xlabel('Temperature (K)')
plt.ylabel('Number of intact molecules')
# plt.legend()
#plt.grid(visible=True)
plt.savefig('..\\python_outputs\\onset\\onset-'+str(random.randint(1000000,9999999)), dpi=1000,bbox_inches='tight')
#-------------------