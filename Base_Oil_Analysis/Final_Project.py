# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Dec 1 09:05:16 2023
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
import random
from sklearn.cluster import KMeans
import pickle
import os

#%% parse bondfile
def parsebondfile(bondfilepath, cutoff=0.3,**kwargs):
    print(f"bond order cutoff {cutoff} is used")
    #10 times faster than version 1
    
    #-------keyward arguments-------------
    bo              = kwargs.get('bo',False) # to get bondorders
    mtypes          = kwargs.get('mtypes',False) # to get molecule ids
    charge          = kwargs.get('charge',False) # get charge of atoms
    abo             = kwargs.get('abo',False)
    nlp             = kwargs.get('nlp',False) # get the number of lone pair
    mols            = kwargs.get('mols',False) # get mtypes wise molecues
    ALL             = kwargs.get('ALL',False) # get everything
    pkl             = kwargs.get('pkl',None) # Fatser the process (directory)
    firstStep_only = kwargs.get('firstStep_only',False) # first step only
    #-------------------------------------
    if firstStep_only:
        pkl = None
        print('If firstStep_only is True, pkl automatically set to None')
        
    if pkl is not None:
        ALL = True
        print('If pkl is not None, ALL automatically set to True')
    #-------------------------------------
    
    
    def parse():
        bonddata = {}
        neighbours = {}
        atomtypes  = {}
        # bondorders[step][atom_id1][atom_id2]=bo between atom1 and atom2
        if bo or ALL: bondorders    = {}
        if mtypes or ALL: molecule_types    = {}
        if charge or ALL: charges = {}
        if nlp or ALL: NLP = {}
        if abo or ALL: ABO = {}
        if mols or ALL: molecules = {}
        if firstStep_only:
            fs_flag = False
        
        with open(bondfilepath) as bf:
            prev_natoms = 0
            natoms_flag = False
            warning_flag = False
            first_warning_ignore = True
            atom_counter = 0
            
            for line in bf:
                splitted = line.split()
                if line.find('Timestep')!=-1:
                    step = int(splitted[-1])
                    
                    if firstStep_only:
                        if fs_flag:
                            break
                        fs_flag = True
                        
                    if bo or ALL: bondorders[step]  = {}
                    if charge or ALL: charges[step] = {}
                    if nlp or ALL: NLP[step] = {}
                    if abo or ALL: ABO[step] = {}
                    
                    neighbours[step]={}
                    atomtypes = {}
                    
                if line.find('Number of particles')!=-1:
                    current_natoms = int(splitted[-1])
                    
                    if atom_counter!=current_natoms:
                        if first_warning_ignore:
                            first_warning_ignore=False
                        else:
                            warning_flag=True
                        
                    if natoms_flag and current_natoms!=prev_natoms:
                        print('User warning from get_neighbours function: Lost atom warning!!')
                    atom_counter = 0   
                    prev_natoms  = current_natoms #new change
                    
                if splitted != [] and splitted[0].isnumeric():
                    atom_counter       +=1
                    parent             = int(splitted[0])
                    parent_type        = int(splitted[1])
                    number_of_children = int(splitted[2])
                    children           = list(map(int,splitted[3:3+number_of_children]))
                    
                    # bond orders
                    bo_id1 = 4+number_of_children
                    bo_idn = bo_id1+number_of_children
                    bo_values = list(map(float,splitted[bo_id1:bo_idn]))
                    
                    ## cutoffs
                    updated_children = []
                    for child,bond_order in zip(children,bo_values):
                        if bond_order>=cutoff:
                            updated_children.append(child)
                    
                    ## assigning values
                    neighbours[step][parent] = updated_children
                    atomtypes[parent]        = parent_type
                    
                    if mtypes or ALL:
                        molecule_types[parent] = int(splitted[3+number_of_children])
                    if bo or ALL:
                        bondorders[step][parent]={}
                        for i,child in enumerate(children):
                            bondorders[step][parent][child] = bo_values[i]
                        
                    if charge or ALL:
                        charges[step][parent] = float(splitted[-1])
                    if nlp or ALL:
                        NLP[step][parent] = float(splitted[-2])
                    if abo or ALL:
                        ABO[step][parent] = float(splitted[-3])
                    if mols or ALL:
                        mol_type = int(splitted[3+number_of_children])
                        if mol_type in molecules.keys():
                            molecules[mol_type].add(parent)
                        else:
                            molecules[mol_type] = set()
            
            if warning_flag:
                print('User warning from get_neighbours function: Repeated atom information detected in some timesteps!')
                
            bonddata['neighbours'] = neighbours
            bonddata['atypes'] = atomtypes
            if bo or ALL: bonddata['bondorders']  = bondorders
            if mtypes or ALL: bonddata['mtypes']  = molecule_types
            if charge or ALL: bonddata['charge']  = charges
            if nlp or ALL: bonddata['nlp'] = NLP
            if abo or ALL: bonddata['abo'] = ABO
            if mols or ALL: bonddata['molecules'] = molecules
            
            if firstStep_only:
                for key, data in bonddata.items():
                    if key in ['neighbours','bondorders']:
                        for step, info in data.items():
                            bonddata[key]=info
                bonddata['fstep'] = step
            
            return bonddata
    
    if pkl is not None:
        if pkl == 'yes':
            path   = path = bondfilepath[:bondfilepath.rfind('.')]+'.pickle'
        else:
            path = pkl
            
        if os.path.exists(path):
            # Load data from pickle
            print('Loading data from pickle...')
            with open(path, 'rb') as pf:
                result = pickle.load(pf)
        else:
            # Execute the provided function to get the data
            print('No pickle file exist. Data is dumping for further use.')
            result = parse()
            with open(path, 'wb') as pf:
                pickle.dump(result, pf)
    
    else:
        result = parse()
        
    return result

#%% bond order evaluation
def bondorder_evolution(bondorders,bondlist,**kwargs):
    if 'ts' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Timestep not specified by the user. The default value of 0.25 fs has been set')
    ## Getting kwargs ##
    ts           = kwargs.get('ts',0.25)
    title        = kwargs.get('title',"")
    savedir      = kwargs.get('savedir', '')
    skipts       = kwargs.get('skipts',None) # skip timestep
    sort         = kwargs.get('sort',None) #ascending,descending or index(ts)
    figsize      = kwargs.get('figsize',(6,10))
    fontsize     = kwargs.get('fontsize',10)
    plot         = kwargs.get('plot','yes')
    ps2temp      = kwargs.get('ps2temp',None) # (it,ramp)
                        # it   = initial temp
                        # ramp = ramping rate
    
    steps = list(bondorders.keys())
    bo_evolution = {}
    ps = np.array(bfp.step2picosecond(steps, ts))
    xlabel = ''#'Time (ps)'
    
    if ps2temp:
        xlabel = 'Temrature (K)'
        it,ramp = ps2temp
        ps = it + ps*ramp

    for index, atompair in enumerate(bondlist):
        atom1,atom2 = atompair
        bo_evolution[index+1]={}
        for step,p in zip(steps,ps):
            bo = bondorders[step][atom1].get(atom2,0)
            bo_evolution[index+1][p]=bo
            
    #############-sort-##################       
    df = pd.DataFrame(bo_evolution)
    df = df.transpose()    
    if sort is not None:
        if sort == 'ascending':
            df.loc[:,'sum'] = df.sum(axis=1)
            df = df.sort_values(by='sum')
            df = df.drop(['sum'],axis=1)      
        elif sort == 'descending':
            df.loc[:,'sum'] = df.sum(axis=1)
            df = df.sort_values(by='sum',ascending=False)
            df = df.drop(['sum'],axis=1)
        elif type(sort) in [float,int]:
            df = df.sort_values(by=sort)
        else:
            print('User Error from bondorder_evolution: sort key not found')
            sys.exit()
    #####################################
    df.index = range(1,len(bondlist)+1)
    #############-skipts-##################
    if skipts is not None:
        df = df.loc[:,skipts:]
        df.columns = df.columns-skipts
    ######################################    
    return df

#%% K-means clustering
baseoil = "PAO4"
colors  = {'PAO4':'blue', 'Squalane':'green'}
directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\Sim-1'.format(baseoil,baseoil)
filename  = '\\bonds.reaxc'
bondfile  = directory+filename
bonddata  = bfp.parsebondfile(bondfile, bo=True, charge=True)

neighbours = bonddata['neighbours']
atypes     = bonddata['atypes']
charges    = bonddata['charge']
bondorders = bonddata['bondorders']
firststepneigh = neighbours[list(neighbours.keys())[0]]

## filtering the data
def get_atomtypes(xxx):
    nd_atypes = []
    for x in xxx:
        nd_atypes.append(atypes[x])
    return np.array(nd_atypes)

for step, atom_charge in charges.items():
    temp_atom   = np.array(list(atom_charge.keys()))
    temp_charge = np.array(list(atom_charge.values()))
    mask   = get_atomtypes(temp_atom) == 2
    atom   = temp_atom[mask]
    charge = temp_charge[mask]
    break # only one step

atom_ids = atom
charges = charge

# Reshape charges for KMeans
charges_reshaped = charges.reshape(-1, 1)

# Apply KMeans clustering
n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(charges_reshaped)

# Custom label based on carbon degree
cluster_labels = kmeans.labels_
labels_meanCharge = []
for i in range(n_clusters):
    labels_meanCharge.append(np.mean(charges[cluster_labels == i]))

mean_charge  =  np.array([])
for lab in cluster_labels:
    mean_charge = np.append(mean_charge, labels_meanCharge[lab])

deg2charge = np.array([-0.30,-0.20,-0.05]) # deg 1, 2, 3
label = np.array([],dtype=int)
for each_charge in mean_charge:
    error    = np.abs(deg2charge-each_charge)
    minErrID = error.argmin()
    label    = np.append(label,minErrID)

atom2label = {}
for lab,atom in zip(label,atom_ids):
    atom2label[atom]=lab
####################-Bonds-##############################
bonds = {}
for parent, children in firststepneigh.items():
    if atypes[parent] != 2:
        continue
    parent_label = atom2label[parent]
    for child in children:
        if atypes[child] != 2:
            continue
        child_label = atom2label[child]
        bond_label  = "".join(sorted(str(parent_label+1)+str(child_label+1)))
        
        bond = tuple(sorted((parent,child)))
        if bond_label in bonds:
            bonds[bond_label].add(bond)
        else:
            bonds[bond_label] = set()


#%% Count the bond and plot
## User input
plot_keys     = sorted(bonds.keys())
skipts        = (1600-300)/4
timestep      = 0.25
steps         = np.array(list(bondorders.keys()))
cbar_min      = 0.0
cbar_max      = 2.0
fontsize      = 15 # All fonts and labels
ticksize      = 10 # Only x-ticks and cbar ticks

# each plot use figsize=(2.25,4)
n_plot = len(plot_keys)
fig, ax = plt.subplots(1,len(plot_keys),figsize=(4.25*n_plot, 4))

stat_reactivity = {}
for i,key in enumerate(plot_keys):
    b = bonds[key]
    df = bfp.bondorder_evolution(bondorders, b, ts=timestep,skipts=skipts,plot='no')
    df_bool = (df<1.1) & (df>0.45)
    df_count= df_bool.sum(axis=0)
    df_norm = (df_count-df_count.min())/(df_count.max()-df_count.min())
    df_norm = df_count
    
    # Fitting a linear function
    x = df_norm.index.values
    y = df_norm.values
    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit)
    
    multiplier = 100
    power      = int(np.log10(multiplier))
    stat_reactivity[key] = abs(multiplier*fit[0]) #taking the slope
    print(fit[0]*multiplier)
    
    # Plotting the data and the fit
    ax[i].scatter(x, y, s=5,label='simulation',color=colors[baseoil])
    ax[i].plot(x, fit_fn(x), color='red',label='fit')
    ax[i].set_title('{}'.format(key))
    ax[i].legend()
    ax[i].set_ylim(0,500)
    ax[i].set_xlim(0,1000)
    # ax[i].set_yticks(range(120,281,40))

# Common X Label
fig.text(0.5, -0.04, 'Time (ps)', ha='center', va='center',fontsize=fontsize)
fig.text(0.08, 0.5, 'Bond count', ha='center', va='center',fontsize=fontsize,rotation=90)

title = directory[directory.find('Base'):]+'\n\n'
fig.suptitle(title, y=1.04, fontsize=10)
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\bondplot'

plt.savefig(savedir+'\\bo_count_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
plt.show()
plt.bar(stat_reactivity.keys(),stat_reactivity.values(),color=colors[baseoil])
plt.xlabel('Types of bonds',fontsize=fontsize-2)
plt.ylabel(f'Statistical bond dissociation rate ($x 10^{power}$)',fontsize=fontsize-4)
plt.title(title,fontsize=10)
plt.ylim(0,12.5)
plt.savefig(savedir+'\\bond_transition_rate_'+str(random.randint(100000, 999999)),dpi=500,bbox_inches='tight')
