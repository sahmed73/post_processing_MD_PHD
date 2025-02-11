# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 23:56:22 2023 (Version 1.0)
Version 2.0 released on June 1, 2023

@author: shihab

This is my own LAMMPS post-processing library
"""

import networkx as nx
from collections import Counter
import re
import sys
import math
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import random
import os
import pickle
import json
from functools import wraps
import numpy as np
from scipy.optimize import curve_fit

# =============================================================================
#   Dacorator functions:
#     1.  function_runtime
#     2.  loadpickle_or_execute
#     3.  loadjson_or_execute (barely used)

#   Module functions:
#     1.  get_neighbours
#     2.  merge_bondfile
#     3.  get_molecules
#     4.  get_molecular_formula      
#     5.  stepwise_species_count (v: get_SpeciesCountAtEveryTimestep)
#     6.  expression_selector
#     7.  pathway_tracker
#     8.  step2picosecond
#     9.  get_nearestindex
#     10. sort_molecular_formula (v: sort_molecular_formula_OLD)
#     11. make_molecular_formula_latex
#     12. compute_molecular_weight
#     13. plot_species_heatmap
#     14. get_speciesVStemp
#     15. get_onset
#     16. bondorder_evolution
# =============================================================================
def function_runtime(f):
    @wraps(f)
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = f(*args,**kwargs)
        stop_time = time.time()
        runtime = stop_time-start_time
        minute = int(runtime/60)
        second = runtime%60
        if minute==0:
            msg = "Execution time of {}: {:0.1f} sec".format(f.__name__,second)
        else:
            msg = "Execution time of {}: {} min {:0.1f} sec".format(f.__name__,minute,second)
        print(msg)
        return result
    return wrapper

def loadpickle_or_execute(pickle_path,function, *args, **kwargs):
    if os.path.exists(pickle_path):
        # Load data from pickle
        start = time.time()
        print('Loading data from pickle instead of {}...'.format(function.__name__))
        with open(pickle_path, 'rb') as pf:
            data = pickle.load(pf)
            stop = time.time()
            _min = int((stop-start)/60)
            _sec = (stop-start)%60
            if _min==0:
                print('Loaded Succesfully!! Loading time: {:0.1f} sec'.format(_sec))
            else:
                print('Loaded Succesfully!! Loading time: {} min {:0.1f} sec'.format(_min,_sec))
            print('-'*60)
    else:
        # Execute the provided function to get the data
        print('Loading data from {}...'.format(function.__name__))
        data = function(*args,**kwargs)
        print('Loaded Succesfully!!!')
        print('-'*60)
        # Save data to pickle for future use
        with open(pickle_path, 'wb') as pf:
            pickle.dump(data, pf)
    return data

def loadjson_or_execute(pickle_path,function, *args, **kwargs):
    if os.path.exists(pickle_path):
        # Load data from pickle
        print('Loading data from JSON instead of {}...'.format(function.__name__))
        print('Loaded Succesfully!!!')
        print('-'*60)
        with open(pickle_path, 'r') as pf:
            data = json.load(pf)
    else:
        # Execute the provided function to get the data
        print('Loading data from {}...'.format(function.__name__))
        data = function(*args,**kwargs)
        print('Loaded Succesfully!!!')
        print('-'*60)
        # Save data to pickle for future use
        with open(pickle_path, 'w') as pf:
            json.dump(data, pf)
    return data

# ---------------------List of Neighbours-----------------------------------
@function_runtime
def get_neighbours(bondfile,**kwargs):
    #10 times faster than version 1
    
    #-------keyward arguments-------------
    bo     = kwargs.get('bo',False) # to get bondorders
    mtypes = kwargs.get('mtypes',False) # to get molecule ids
    charge = kwargs.get('charge',False) # get charge of atoms
    nlp    = kwargs.get('nlp',False) # get the number of lone pair
    #-------------------------------------
    
    neighbours = {}
    atomtypes  = {}
    if bo: bondorders    = {}
    if mtypes: molecule_types    = {}
    if charge: charge = {}
    if nlp: nlp = {}
    
    with open(bondfile) as bf:
        prev_natoms = 0
        natoms_flag = False
        warning_flag = False
        first_warning_ignore = True
        atom_counter = 0
        
        for line in bf:
            splitted = line.split()
            if line.find('Timestep')!=-1:
                step = int(splitted[-1])
                neighbours[step]={}
                atomtypes = {}
                if bo: bondorders[step]  = {}
                if charge : charge[step] = {}
                if nlp: nlp[step] = {}
                
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
            if splitted != [] and splitted[0].isnumeric():
                atom_counter            +=1
                parent                   = int(splitted[0])
                parent_type              = int(splitted[1])
                number_of_children       = int(splitted[2])
                children                 = list(map(int,splitted[3:3+number_of_children]))
                neighbours[step][parent] = children
                atomtypes[parent]        = parent_type
                
                if mtypes:
                    molecule_types[parent] = int(splitted[3+number_of_children])
                if bo:
                    bo_id = 4+number_of_children
                    bo_values = list(map(float,splitted[bo_id:bo_id+number_of_children]))
                    bondorders[step][parent]={}
                    for i,child in enumerate(children):
                        bondorders[step][parent][child] = bo_values[i]
                if charge:
                    charge[step][parent] = float(splitted[-1])
                if nlp:
                    nlp[step][parent] = float(splitted[-2])
        
        if warning_flag:
            print('User warning from get_neighbours function: Repeated atom information detected in some timesteps!')
            
        result = (neighbours,atomtypes,)
        if bo: result+=(bondorders,)
        if mtypes: result+=(molecule_types,)
        
        return result

# ---------------------List of Neighbours-----------------------------------
@function_runtime
def parsebondfile(bondfilepath,**kwargs):
    #10 times faster than version 1
    
    #-------keyward arguments-------------
    bo     = kwargs.get('bo',False) # to get bondorders
    mtypes = kwargs.get('mtypes',False) # to get molecule ids
    charge = kwargs.get('charge',False) # get charge of atoms
    nlp    = kwargs.get('nlp',False) # get the number of lone pair
    mols   = kwargs.get('nlp',False) # get mtypes wise molecues
    ALL    = kwargs.get('ALL',False) # get everything
    pkl    = kwargs.get('pkl',None) # Fatser the process (directory)
    fso    = kwargs.get('fso',False) # first step only
    freq   = kwargs.get('freq',False) # frequency
    #-------------------------------------
    if fso:
        pkl = None
        print('If fso is True, pkl automatically set to None')
        
    if pkl is not None:
        ALL = True
        print('If pkl is not None, ALL automatically set to True')
    #-------------------------------------
    
    def parse():
        bonddata = {}
        neighbours = {}
        atomtypes  = {}
        if bo or ALL: bondorders    = {}
        if mtypes or ALL: molecule_types    = {}
        if charge or ALL: charges = {}
        if nlp or ALL: NLP = {}
        if mols or ALL: molecules = {}
        if fso:
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
                    
                    if fso:
                        if fs_flag:
                            break
                        fs_flag = True
                        
                    if bo or ALL: bondorders[step]  = {}
                    if charge or ALL: charges[step] = {}
                    if nlp or ALL: NLP[step] = {}
                    
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
                    
                if splitted != [] and splitted[0].isnumeric():
                    atom_counter            +=1
                    parent                   = int(splitted[0])
                    parent_type              = int(splitted[1])
                    number_of_children       = int(splitted[2])
                    children                 = list(map(int,splitted[3:3+number_of_children]))
                    neighbours[step][parent] = children
                    atomtypes[parent]        = parent_type
                    
                    if mtypes or ALL:
                        molecule_types[parent] = int(splitted[3+number_of_children])
                    if bo or ALL:
                        bo_id = 4+number_of_children
                        bo_values = list(map(float,splitted[bo_id:bo_id+number_of_children]))
                        bondorders[step][parent]={}
                        for i,child in enumerate(children):
                            bondorders[step][parent][child] = bo_values[i]
                    if charge or ALL:
                        charges[step][parent] = float(splitted[-1])
                    if nlp or ALL:
                        NLP[step][parent] = float(splitted[-2])
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
            if bo or ALL: bonddata['bondorders'] = bondorders
            if mtypes or ALL: bonddata['mtypes'] = molecule_types
            if charge or ALL: bonddata['charge'] = charge
            if nlp or ALL: bonddata['nlp'] = nlp
            if mols or ALL: bonddata['molecules'] = molecules
            
            if fso:
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

@function_runtime
def merge_bondfiles(*args,**kwargs):
    # input: bondfilepath_1, start_ts_1, end_ts_1,
    #       bondfilepath_2, start_ts_2, end_ts_2, ...., **kwargs
    
    #-------keyward arguments-------------
    file     = kwargs.get('file','merge_bonds.reaxc')
    # nevery   = kwargs.get('nevery',None)
    #-------------------------------------
    
    #------arguments----------------------
    paths, start_ts, end_ts = [], [], []
    for i,arg in enumerate(args):
        if i%3==0:
            paths.append(arg)
        elif i%3==1:
            start_ts.append(arg)
        else:
            end_ts.append(arg)
    #-------------------------------------
    with open(file,'w') as mf:
        for i, bondfilepath in enumerate(paths):
            what_todo = 'do nothing' # flag for start writing
            with open(bondfilepath,'r') as bf:
                sts = start_ts[i]
                ets = end_ts[i]
                
                for line in bf:
                    if line.find('Timestep')!=-1:
                        cts = int(line.split()[-1]) # current timestep
                        if cts == sts:
                            what_todo = 'start writing'
                        if cts == ets:
                            what_todo = 'stop writing & break'
                    
                    if what_todo == 'start writing':
                        mf.write(line)
                    if what_todo == 'stop writing & break':
                        break
                


def merge_bondfiles_OLD(*args,**kwargs):
    
    #-------keyward arguments-------------
    file = kwargs.get('file','marge_bonds.reaxc')
    #-------------------------------------
    
    startline = []
    stopline  = []
    for index,arg in enumerate(args):
        skip_first = True
        with open(arg,'r') as bf:
            for line in bf:
                if index == 0 and line.find('Timestep')!=-1:
                    startline.append(line)
                    break
                if line.find('Timestep')!=-1:
                    if skip_first:
                        skip_first = False
                        continue
                    startline.append(line)
                    break
                
    for i in range(len(startline)-1):
        stopline.append(startline[i+1])
    stopline.append(-1) # till end
    
    with open(file,'w') as mf:
        for index,arg in enumerate(args):
            start_writing = False
            with open(arg,'r') as bf:
                for line in bf:
                    if line == stopline[index] and stopline[index]!=-1:
                        break
                    if not line.endswith('\n'):
                        line+='\n'
                    if start_writing:
                        mf.write(line)
                    if line == startline[index]:
                        if line.find('Timestep')!=-1:
                            mf.write(line)
                            start_writing = True
    print('Merged bond file created!')

#--------------Is equal two text file -----------------------------------------                           
def text_isequal(textfile1,textfile2,strip=False):
    result = True
    with open(textfile1,'r') as tf1, open(textfile2,'r') as tf2:
        for i,lines in enumerate(zip(tf1,tf2)):
            line1,line2=lines
            if strip:
                line1 = line1.strip()
                line2 = line2.strip()
            if line1!=line2:
                result = False
                print('Mismatched:', 'line',i+1)
    return result
    
    
# ---------------------get molecule from list neighbour------------------------------
def get_molecules(neigh):
    '''
    Parameters
    ----------
    neigh : Dictionary
        DESCRIPTION.
        neighbours of a single timestep or neighbours[step]

    Returns
    -------
    molecules : Set
        DESCRIPTION.
    '''
    graph = nx.Graph(neigh)
    molecules = nx.connected_components(graph)    
    return molecules

# -----------------------get chemical formula-------------------------------------
def get_molecular_formula(molecule,atomtypes,atomsymbols):
    '''
    Parameters
    ----------
    molecule : any itterable
        DESCRIPTION.
    atomtypes : TYPE
        DESCRIPTION.
    atomsymbols : TYPE, optional
        DESCRIPTION. The default is 'HCO'.

    Returns
    -------
    formula : TYPE
        DESCRIPTION.

    '''
    species = map(lambda x: atomtypes[x],molecule)
    counter = [0]*len(atomsymbols)
    for s in species:
        counter[s-1]+=1
    
    formula = ''
    for i in range(len(atomsymbols)):
        if counter[i]!=0:
            if counter[i]!=1:
                formula += atomsymbols[i]+str(counter[i])
            else:
                formula += atomsymbols[i]
    return formula

# ------------Get All Speceis Count At every Timestep--------
@function_runtime
def get_SpeciesCountAtEveryTimestep(neighbours,atomtypes,atomsymbols,exception=[],step2ps=False,step2psargs=[],minps=-math.inf,maxps=math.inf):
    '''
    Parameters
    ----------
    neighbours : Dictionary
        DESCRIPTION.
    atomtypes : Dictionary
        DESCRIPTION.
    atomsymbols : any iterable
        DESCRIPTION.
    exception : list (any iterable), optional
        DESCRIPTION. The default is [].
    step2ps : bool, optional
        DESCRIPTION. The default is False.
    step2psargs : list (any iterable), optional
        DESCRIPTION. The default is [].
    minps : TYPE, optional
        DESCRIPTION. The default is -math.inf.
    maxps : TYPE, optional
        DESCRIPTION. The default is math.inf.

    Returns
    -------
    stepwise_species_count : Dictionary
        DESCRIPTION.

    '''
    #local function
    def get_species(molecules,atomtypes,atomsymbols): #take the molecule file from get_molecute()
        #atomlist=['H','C','O'], same order as TYPE
        allspecies = []
        for molecule in molecules:
            species = get_molecular_formula(molecule, atomtypes,atomsymbols)
            allspecies.append(species)
        allspecies_count = Counter(allspecies)
        
        return allspecies_count

    ## if step2ps=True --> Convert timesteps to piccoseconds ##
    _steps  = list(neighbours.keys())
    s2ps = _steps.copy()
    s2ps_dict = {}
    
    if step2ps:
        if not step2psargs:
            sys.exit('Value error: \'step2psargs\' is empty')
        del s2ps
        s2ps = step2picosecond(_steps, *step2psargs)
    
    s2ps_dict = dict(zip(_steps,s2ps))
    ##-------------------------------------------------------##
    
    stepwise_species_count = {}
    for step,neigh in neighbours.items():
        molecules = get_molecules(neigh)
        species_count = get_species(molecules, atomtypes,atomsymbols)
        
        ## Delete species from the Exception list ##
        for ex in exception: del species_count[ex]
        ##------------------------------------------##
        if minps<=s2ps_dict[step]<maxps:
            stepwise_species_count[s2ps_dict[step]]=species_count
                
    return stepwise_species_count     

@function_runtime
def stepwise_species_count(neighbours,atomtypes,atomsymbols,step2ps=None):
    # if requires ps instead of steps set step2ps = timestep
    _swsc = {} # stepwise specis count _swsc[step][species]=count
    if step2ps:
        steps = list(neighbours.keys())
        ps    = step2picosecond(steps, step2ps)
    for i,items in enumerate(neighbours.items()):
        step,neigh = items
        molecules = get_molecules(neigh)
        _sc = {} # species count
        for molecule in molecules:
            species = get_molecular_formula(molecule, atomtypes,atomsymbols)
            if species in _sc.keys():
                _sc[species]+=1
            else:
                _sc[species] = 1
                
        if step2ps:
            _swsc[ps[i]]=_sc
        else:
            _swsc[step]=_sc
    return _swsc    

@function_runtime
def stepwise_species_count_v2(atomtypes,atomsymbols,*args,**kwargs):
    # if requires ps instead of steps set step2ps = timestep
    step2ps = kwargs.get('step2ps',None)
    kind    = kwargs.get('kind','sum')
    
    def inner(neighbours):
        # if requires ps instead of steps set step2ps = timestep
        _swsc = {} # stepwise specis count _swsc[step][species]=count
        if step2ps:
            steps = list(neighbours.keys())
            ps    = step2picosecond(steps, step2ps)
        for i,items in enumerate(neighbours.items()):
            step,neigh = items
            molecules = get_molecules(neigh)
            _sc = {} # species count
            for molecule in molecules:
                species = get_molecular_formula(molecule, atomtypes,atomsymbols)
                if species in _sc.keys():
                    _sc[species]+=1
                else:
                    _sc[species] = 1
                    
            if step2ps:
                _swsc[ps[i]]=_sc
            else:
                _swsc[step]=_sc
        return _swsc
    
    result = []
    for arg in args:
        result.append(inner(arg))
    
    output = {}
    for out in result:
        for step, species_count in out.items():
            if step not in output.keys():
                output[step]={}
            for species,count in species_count.items():
                if species not in output[step].keys():
                    output[step][species] = count
                else:
                    output[step][species] +=count
    if kind=='sum':
        pass
    elif kind=='mean':
        n = len(args)
        print(n)
        for step, species_count in output.items():
            for species,count in species_count.items():
                    cc = int(count/n)
                    output[step][species] = cc+1 if cc!=count/n else cc
    else:
        print('User Error from stepwise_species_count: kind key error')
    return output
        

# ---------------------Expression Selector----------------------------------------
def expression_selector(neighbours,atomtypes,file='expression_selector.txt',atomsybols='HCO',heading=None):
    steps = neighbours.keys()
    with open(file, 'w') as f:
        if heading:
            f.write(heading+'\n'+'-------------------------------'+'\n')
        for step in steps:
            molecules = get_molecules(neighbours[step])
            f.write('Timestep='+str(step)+'\n')
            
            for molecule in molecules:
                formula = get_molecular_formula(molecule, atomtypes,atomsybols)
                f.write('Timestep='+str(step)+'\t')
                f.write('Molecule: '+formula+'\n')
                f.write('Molecule: '+str(molecule)+'\n')
                for atom in molecule:
                    f.write('ParticleIdentifier=='+str(atom)+'|| ')
                f.write('\n\n')
                
# ---------------Pathway Tracker Of a singler molecule---------------------------
def pathway_tracker(seek_molecule,neighbours,atomtypes,file='pathway_tracker.txt',atomsybols='HCO'):
    seek_molecule_formula = get_molecular_formula(seek_molecule, atomtypes)
    steps = neighbours.keys()
    file = '{}'.format(seek_molecule)+seek_molecule_formula+'_'+file
    
    with open(file, 'w') as f:
        f.write('This is a pathway tracker for {}'.format(seek_molecule_formula)+'\n')
        f.write('{}'.format(seek_molecule))
        f.write('-----------------------------------------------------------------\n\n')
        pathway = []
        for step in steps:
            molecules = get_molecules(neighbours[step])
            for molecule in molecules:
                formula = get_molecular_formula(molecule,atomtypes,'HCO')
                intersection = set(seek_molecule) & set(molecule)
                if intersection and formula not in pathway:
                    pathway.append(formula)
                    f.write('Timestep='+str(step)+'\t')
                    f.write('Molecule: '+formula+'\n')
                    for atom in molecule:
                        f.write('ParticleIdentifier=='+str(atom)+'|| ')
                        
                    f.write('\n{}\n\n'.format(molecule))
    return pathway

# --------------------------------------------------------------
def step2picosecond(steps,*args):
    timestep =[]
    if len(args)==1:
        tstep1, = args
        current_time = steps[0]*tstep1/1000
        previous_step = steps[0]
        
        for current_step in steps:
            current_time += (current_step-previous_step)*tstep1/1000  
            previous_step = current_step
            timestep.append(current_time)
        return timestep
    
    elif len(args)==3:
        tstep1,tstep2,limit = args
        current_time = steps[0]*tstep1/1000
        previous_step = steps[0]
        
        for current_step in steps:
            if current_step<=limit: 
                current_time += (current_step-previous_step)*tstep1/1000
            else:        
                current_time += (current_step-previous_step)*tstep2/1000
            previous_step = current_step
            timestep.append(current_time)
        return timestep
    else:
        sys.exit('Error: The length of the argument in the step2picosecond() function is not specified correctly.\n')


# -----------getting index of specific timestep---------------
def get_nearestindex(_timestep,*args,**kwargs):
    if len(args)==1:
        _time = args[0]
    key = kwargs.get('key',None) # posible value: 'first','last'
    if key=='first':
        index = 0
    elif key=='last':
        index = len(_timestep)-1
    elif key==None:
        _temp = [abs(_value-_time) for _value in _timestep]
        _mini = min(_temp)
        index = _temp.index(_mini)
    else:
        sys.exit('User Error from get_nearestindex: Key not found!')
    return index


# -----------sort chemical formula in this order: C,H,O---------------
def sort_molecular_formula(molecular_formula,order=['C','H','O']):
    def do_sort(species):
        item = re.findall('[A-Z][a-z]?|\d+|.', species)+['']
        elements = []
        for i in range(len(item)-1):
            if item[i].isalpha() and item[i+1].isnumeric():
                elements.append((item[i],item[i+1]))
            elif item[i].isalpha() and not item[i+1].isnumeric():
                elements.append((item[i],''))
                
        symbols = [x for x,_ in elements]
        if not set(symbols).issubset(set(order)):
            sys.exit('User Error: Some elements does not in the order list! Please specify order list accordingly.')
        
        for i in range(len(elements)):
            elements[i]+=(order.index(elements[i][0]),)
        elements.sort(key=lambda x:x[2])
        sorted_chem_formula = ''
        for element in elements:
            sorted_chem_formula+=element[0]+element[1]
        return sorted_chem_formula
          
    
    if type(molecular_formula)!=str:
        result = []
        for species in molecular_formula:
            result.append(do_sort(species))
    else:
        result = do_sort(molecular_formula) 
    return result
   
# ------------------make_latex-----------------------------------------
def make_molecular_formula_latex(molecular_formula,sort=False):
    
    if sort:
        formula = list(sort_molecular_formula(molecular_formula))
    else:
        formula = list(molecular_formula)
    
    
    latex_formula = ['']*len(formula)
    for i in range(len(formula)):
        formula[i]+='$'
        for j in range(len(formula[i])-1):
            if formula[i][j].isalpha() and formula[i][j+1].isnumeric():
                latex_formula[i]+=formula[i][j] + '_{'
            elif formula[i][j].isnumeric() and not formula[i][j+1].isnumeric():
                latex_formula[i]+=formula[i][j] + '}'
            else:
                latex_formula[i]+=formula[i][j]
        #if formula[i][-2].isnumeric(): latex_formula[i]+='}'
        latex_formula[i] = '$'+latex_formula[i]+'$'
    
    return latex_formula

# --------------------------------------------------------------
def compute_molecular_weight(species,exclude=[]):
    atomic_info    = atomic_weight(key='all')
    item = re.findall('[A-Z][a-z]?|\d+|.', species)
    item.append('q') #fake
    molecular_weight = 0.0
    for i in range(len(item)-1):
        now,nxt = item[i],item[i+1]
        if now in exclude:
            continue
        if now.isalpha() and nxt.isalpha():
            if now in atomic_info:
                molecular_weight+=atomic_info[now]
            else:
                sys.exit('Error: No such element \'{}\' in the compute_molecular_weight() function'.format(now))
        elif now.isalpha() and nxt.isnumeric():
            if now in atomic_info:
                molecular_weight+=atomic_info[now]*float(nxt)
            else:
                sys.exit('Error: No such element \'{}\' in the compute_molecular_weight() function'.format(now))
    return molecular_weight

# ----------------------Heatmap of Species--------------------------------
@function_runtime
def plot_species_heatmap(atomtypes,atomsymbols,*neighbours_list,**kwargs):
    # keyward agruments:
        # nspecies (int)   = number of species to be plot (default value = 20)
        # ts (int)         = timestep (default value = 0.25)
        # savedir (string) = directory to save (default value = '')
        # titile (string)  = title to be shown on the top (default value = '')
        # order (list)     = order of the species. Takes list of species (default = mean abundance)
        # pikle (string)   = load from pickle. It takes pickle_path (default = None)
        # skipts (float)   = skip some data from the first (default = None)
        # topspec (list)   = species that appears on the top (default = None)
        # 
        
    nspecies    = kwargs.get('nspecies', 20)
    ts          = kwargs.get('ts', 0.25)
    savedir     = kwargs.get('savedir', '')
    title       = kwargs.get('title', '')
    order       = kwargs.get('order',None)
    pickle      = kwargs.get('pickle',None)
    skipts      = kwargs.get('skipts',None) # skip timestep
    topspec     = kwargs.get('topspec',None)
    ignor       = kwargs.get('ignor',None)
    kind        = kwargs.get('kind','sum') # kind = sum,mean
    log         = kwargs.get('log',True)
    exclude     = kwargs.get('exclude',None)
    fontsize    = kwargs.get('fontsize', 12)
    figsize     = kwargs.get('figsize',(15,8))
    
    if pickle:
        data = loadpickle_or_execute(pickle, stepwise_species_count_v2, atomtypes, atomsymbols, *neighbours_list, step2ps=ts)
    else:
        data = stepwise_species_count_v2(atomtypes, atomsymbols,*neighbours_list,step2ps=ts,kind=kind)
    df = pd.DataFrame(data).fillna(1)
    
    if order:
        df = df.loc[order,:]
    else:
        df.loc[:,'sum'] = df.sum(axis=1)
        df = df.sort_values(by='sum',ascending=False)
        df = df.drop(['sum'],axis=1)
        
    #############-exclude-##################
    if exclude is not None:
        df = df.drop(index=exclude)
    ######################################
    
    #############-skipts-##################
    if skipts is not None:
        df = df.loc[:,skipts:]
        df.columns = df.columns-skipts
    ######################################
    
    #############-ignor-##################
    if ignor is not None:  # ignor = (main_species,'close')
        if isinstance(ignor, tuple):
            main_species, kind = ignor
            main = compute_molecular_weight(main_species)
            
            H = compute_molecular_weight('H')
            O = compute_molecular_weight('O')
            criteria = []
            
            for a in [0,1,2,3]:
                for b in [0,1,2,3]:
                    if a==0 and b==0:
                        continue
                    cri = a*H+b*O
                    criteria.append(cri)
            
            for chem in df.index:
                current = compute_molecular_weight(chem)            
                for cri in criteria:
                    if abs(abs(main-current)-cri)<=0.01:
                        print(chem,'ignored')
                        df = df.drop(chem)
                        continue
        elif isinstance(ignor, str):
            df = df.drop([ignor],axis=0)
    #####################################
    
    
    df = df.head(nspecies)
    ##############-topspec-###############
    if topspec is not None:
        desired_order = []
        for spec in topspec:
            if spec in df.index:
                desired_order.append(spec)
        for spec in df.index:
            if spec not in desired_order:
                desired_order.append(spec)
        df = df.reindex(index=desired_order)
    ######################################
        
    print('Species ORDER:',list(df.index))
    df.index = make_molecular_formula_latex(df.index,sort=True)
    _, ax1 = plt.subplots(figsize=figsize)
    if log:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                     cbar_kws={'label': 'Number of species'},
                     norm=LogNorm(),xticklabels=1000)
    else:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                     cbar_kws={'label': 'Number of species'},
                     xticklabels=1000)
    ax.set_xlabel('Time (ps)',fontsize=fontsize+1)
    # ax.set_ylabel('Species')
    plt.tick_params(left=False,bottom=False)
    plt.yticks(rotation=0,fontsize=fontsize)
    fig = ax.get_figure()
    title += '-{}\nTop {} species\n'.format(kind,nspecies)
    ax.set_title(title)
    ax.set_xticklabels([0,250,500,750,1000],fontsize=fontsize)
    
    
    # cbar
    ax.figure.axes[-1].yaxis.label.set_size(fontsize+1) # cbar label size
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=fontsize)
    
    saveplot = 'species_heatmap_'+str(random.randint(100000, 999999))
    if savedir: savedir += '\\'+saveplot
    else: savedir += saveplot
    fig.savefig(savedir, dpi=400, bbox_inches='tight')
    
# ----------------------Heatmap of Species from multiple neighbours--------------------------------
@function_runtime
def plot_species_heatmap_v2(neighbours,atomtypes,atomsymbols,**kwargs):
    # keyward agruments:
        # nspecies (int)   = number of species to be plot (default value = 20)
        # ts (int)         = timestep (default value = 0.25)
        # savedir (string) = directory to save (default value = '')
        # titile (string)  = title to be shown on the top (default value = '')
        # order (list)     = order of the species. Takes list of species (default = mean abundance)
        # pikle (string)   = load from pickle. It takes pickle_path (default = None)
        # skipts (float)   = skip some data from the first (default = None)
        # topspec (list)   = species that appears on the top (default = None)
        
    nspecies    = kwargs.get('nspecies', 20)
    ts          = kwargs.get('ts', 0.25)
    savedir     = kwargs.get('savedir', '')
    title       = kwargs.get('title', '')
    order       = kwargs.get('order',None)
    pickle      = kwargs.get('pickle',None)
    skipts        = kwargs.get('skipts',None) # skip timestep
    topspec     = kwargs.get('topspec',None)
    ignor       = kwargs.get('ignor',None)
    kind        = kwargs.get('kind','sum') # kind = 'sum' or 'mean'
    log         = kwargs.get('log',True)
    exclude     = kwargs.get('exclude',None)
    
    if pickle:
        data = loadpickle_or_execute(pickle, stepwise_species_count, neighbours, atomtypes, atomsymbols,step2ps=ts)
    else:
        data = stepwise_species_count(neighbours, atomtypes, atomsymbols,
                                      step2ps=ts)
    ###########-fill with###############
    if log:
        df = pd.DataFrame(data).fillna(1)
    else:
        df = pd.DataFrame(data).fillna(0)
    ####################################
    
    if order:
        df = df.loc[order,:]
    else:
        df.loc[:,'sum'] = df.sum(axis=1)
        df = df.sort_values(by='sum',ascending=False)
        df = df.drop(['sum'],axis=1)
    
    #############-exclude-##################
    if exclude is not None:
        df = df.drop(index=exclude)
    ######################################
    
    #############-skipts-##################
    if skipts is not None:
        df = df.loc[:,skipts:]
        df.columns = df.columns-skipts
    ######################################
    
    #############-ignor-##################
    if ignor is not None:  # ignor = (main_species,'close')
        if isinstance(ignor, tuple):
            main_species, kindd = ignor
            main = compute_molecular_weight(main_species)
            
            H = compute_molecular_weight('H')
            O = compute_molecular_weight('O')
            criteria = []
            
            for a in [0,1,2,3]:
                for b in [0,1,2,3]:
                    if a==0 and b==0:
                        continue
                    cri = a*H+b*O
                    criteria.append(cri)
            
            for chem in df.index:
                current = compute_molecular_weight(chem)            
                for cri in criteria:
                    if abs(abs(main-current)-cri)<=0.01:
                        print(chem,'ignored')
                        df = df.drop(chem)
                        continue
        elif isinstance(ignor, str):
            df = df.drop([ignor],axis=0)
    #####################################
    
    
    df = df.head(nspecies)    
    ##############-topspec-###############
    if topspec is not None:
        desired_order = []
        for spec in topspec:
            if spec in df.index:
                desired_order.append(spec)
        for spec in df.index:
            if spec not in desired_order:
                desired_order.append(spec)
        df = df.reindex(index=desired_order)
    ######################################
        
    print('Species ORDER:',list(df.index))
    df.index = make_molecular_formula_latex(df.index,sort=True)
    #df.loc['$O_{2}$',:] = df.loc['$O_{2}$',:] + 200
    _, ax1 = plt.subplots(figsize=(15,8))
    if log:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                         cbar_kws={'label': 'Number of species'},
                         norm=LogNorm(),xticklabels=1000)
    else:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                         cbar_kws={'label': 'Number of species'},
                         xticklabels=1000)
        
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Species')
    plt.tick_params(left=False,bottom=False)
    plt.yticks(rotation=0)
    fig = ax.get_figure()
    title += '---{}\nTop {} species'.format(kind,nspecies)
    ax.set_title(title)
    plt.show()
    
    saveplot = 'species_heatmap_'+str(random.randint(100000, 999999))
    if savedir: savedir += '\\'+saveplot
    else: savedir += saveplot
    fig.savefig(savedir, dpi=400, bbox_inches='tight')
# ------Get Atomic Mass --------------------------------
def atomic_weight(*args,**kwargs):
    # args   = element string
    # kwargs = {key:'all'}
    key = kwargs.get('key',None)
    MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182,
                      'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,'F': 18.9984032,
                      'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
                      'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
                      'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415,
                      'Cr': 51.9961, 'Mn': 54.938045, 'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934,
                      'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216,
                      'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62,
                      'Y': 88.90585, 'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063,
                      'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411,
                      'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6, 'I': 126.90447,
                      'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547,
                      'Ce': 140.116, 'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36,
                      'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5,
                      'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421,
                      'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49,
                      'Ta': 180.9479, 'W': 183.84, 'Re': 186.207,
                      'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
                      'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833,
                      'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.9824,
                      'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197,
                      'Ra': 226.0254, 'Ac': 227.0278, 'Th': 232.03806,
                      'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482,
                      'Pu': 244.0642, 'Am': 243.0614, 'Cm': 247.0703,
                      'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829,
                      'Fm': 257.0951, 'Md': 258.0951, 'No': 259.1009,
                      'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271,
                      'Bh': 270, 'Hs': 269, 'Mt': 278, 'Ds': 281,
                      'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289,
                      'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,'': 0}

    if len(args)==1:
        element = args[0]
        return MM_of_Elements[element]
    elif key=='all':
        return MM_of_Elements
    else:
        sys.exit('User Error from atomic_weight: Unknown kwargs or number of args')
# ------ONSET Finding----------------------
@function_runtime
def get_speciesVStemp(species,neighbour,atomtypes,atomsymbols,**kwargs):
    ########## Soft Warnings ################
    if 'ramp' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Ramping rate not specified by the user. The default value of 1 K/ps has been set')
    if 'ts' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Timestep not specified by the user. The default value of 0.25 fs has been set')
    if 'it' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Initial temperature not specified by the user. The default value of 0 K has been set')
        
    ######### Getting kwargs #################    
    ramp         = kwargs.get('ramp',1) # temp ramp rate (Default: 1 K/ps)
    it           = kwargs.get('it',0) # initial temp (Default: 0 K)
    ts           = kwargs.get('ts',0.25)
    method       = kwargs.get('method','mf')
                   # calculate onset based on method
                   # mf = molecular formula
                   # mw = molecular weight
    lim        = kwargs.get('lim',0)
                   # if method='mw', it will count species if the mw
                   # within the 'lim'
                   
    #-----Fit Function----------------------
    def fit(x, A, B, C, D):   
        y = C+((D-C)/(1+np.exp((x-A)/B)))
        return y

    def inv_fit(y,A,B,C,D):
        x = A + B*np.log((D-y)/(y-C))
        return x
    #-----------------------------------------
    
    if method=='mf':
        count = {}
        for step,neigh in neighbour.items():
            count[step]=0
            molecules = get_molecules(neigh)
            for molecule in molecules:
                mf = get_molecular_formula(molecule, atomtypes, atomsymbols)
                if mf==species:
                    #print(mf)
                    count[step]+=1
            #print(step,count[step])

    elif method=='mw':
        count = {}
        mw_species = compute_molecular_weight(species)
        for step,neigh in neighbour.items():
            count[step]=0
            molecules = get_molecules(neigh)
            for molecule in molecules:
                mf = get_molecular_formula(molecule, atomtypes, atomsymbols)
                mw = compute_molecular_weight(mf)
                if abs(mw-mw_species)<=lim:
                    count[step]+=1
                    #print(mw,mw_species,end='\t')
    else:
        sys.exit('Method not found!!!')
    
    steps          = list(count.keys())
    sc             = np.array(list(count.values())) # species count
    ps             = np.array(step2picosecond(steps, ts))
    temp           = it + ramp*ps
    return temp,sc

def get_onset(temp,sc,**kwargs):
    imc          = kwargs.get('imc',50) # initial molecular count
    ig           = kwargs.get('ig',[2200,107,0,50]) # initial guess for curve fit
    show_fit     = kwargs.get('show_fit',False)
    
    #-----Fit Function----------------------
    def fit(x, A, B, C, D):   
        y = C+((D-C)/(1+np.exp((x-A)/B)))
        return y

    def inv_fit(y,A,B,C,D):
        if D-y<0:
            print('Warning from get_onset: Negative D value adjusted!!')
            D = 49.98
        x = A + B*np.log((D-y)/(y-C))
        return x
    #-----------------------------------------
    
    p, cov = curve_fit(fit, temp, sc,ig,maxfev=2000)
    if show_fit :
        print('Parameters: ',p)
        print('Covariance',cov)
    sc_fit         = fit(temp,*p)
    onset          = inv_fit(imc-1, *p) # onset temp at which 1 molecule
    
    return onset, sc_fit

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
    ps = np.array(step2picosecond(steps, ts))
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
    
    if plot=='yes':
        _, ax1 = plt.subplots(figsize=figsize)
        hm = sns.heatmap(df,cmap='jet',cbar_kws={ 'label': 'Bond order',},xticklabels=1000,ax=ax1,vmin=0.0,vmax=3.0)
        hm.set_xlabel(xlabel,fontsize=fontsize+2)
        # hm.set_ylabel('Bonds ',fontsize=fontsize)
        plt.yticks(rotation=0)
        
        #############-Y Tick Labels-###########
        hm.set_yticklabels([])
        
        ## set xticks location and labels ####
        # hm.set_xticks([795*x for x in range(5)])
        print([int(x) for x in hm.get_xticks()])
        hm.set_xticklabels([0,250,500,750,1000])
        plt.xticks(rotation=90,fontsize=fontsize+2)
        ###########
        
        hm.set_title(title,fontsize=fontsize-5)
        # Set the font size of the y-label
        hm.set_yticklabels(hm.get_yticklabels(), fontsize=fontsize)
    
    
        ######## colorbar settings #########
        cbar = hm.collections[0].colorbar
        cbar.ax.tick_params(labelsize=fontsize+2)
        hm.figure.axes[-1].yaxis.label.set_size(15+4)
        ####################################
        
        if savedir:
            fig = hm.get_figure()
            savedir = savedir+'\\bo_evolution_'+str(random.randint(100000, 999999))
            fig.savefig(savedir, dpi=400, bbox_inches='tight')
    
    return df

@function_runtime
def get_nbondsVStime(bondorders, atomtypes, bondlist,**kwargs):
    # plot the number of bonds vs time
    ## Getting kwargs and setting up default values ##
    cutoff = kwargs.get('cutoff',1)
    tol    = kwargs.get('tol',0.20)
    ts     = kwargs.get('ts',0.25) # timestep
    plot   = kwargs.get('plot',None) # values: 'yes'
    skipts = kwargs.get('skipts',None) # skip timestep
    cumu   = kwargs.get('cumu',False) 
    
    ## Getting args
    
    nbonds = []
    cumu_nbonds = []
    cumu_count = 0
    steps  = []
    for step,bo in bondorders.items():
        steps.append(step)
        bondcount = 0
        for bond in bondlist:
            atom1, atom2 = bond
            try:
                bovalue = bo[atom1][atom2]
            except:
                bovalue = 0
                
            # if bovalue!=0: print(step,bovalue,end=' | ')
            if abs(bovalue-cutoff)<=tol:
                bondcount+=1
                cumu_count+=1
        nbonds.append(bondcount)
        cumu_nbonds.append(cumu_count)
        
    ps = step2picosecond(steps, ts)
    
    if skipts is not None:
        startid = ps.index(skipts)
        ps      = np.array(ps[startid+1:])
        nbonds  = np.array(nbonds[startid+1:])
        if cumu:
            cumu_nbonds  = np.array(cumu_nbonds[startid+1:])
        ps = ps - skipts
        
    # ploting
    if plot == 'yes':
        plt.plot(ps,nbonds,color='black')
        plt.xlabel('Time (ps)')
        plt.ylabel('Number of bonds')
    
    if cumu:
        return ps, cumu_nbonds
    else:
        return ps, nbonds
        
    

@function_runtime
def species_to_molecule(bonddata,atomsymbols, species,**kwargs): # species-->molecule
    # speceis is chemcical formula of a molecule
    # molecule here is the list of atom ids that refer the molecule
    
    ###########-get kwargs-###############
    source    = kwargs.get('source',False)
    dump      = kwargs.get('dump',None)
    shortinfo = kwargs.get('shortinfo',False)
    ######################################
    def atomids2expression(atoms): # atom ids to ovito expression selection
        atoms = list(atoms)
        expression = ''
        for atom in atoms:
            if atoms.index(atom)==len(atoms)-1:
                expression+='ParticleIdentifier=='+str(atom)
            else:
                expression+='ParticleIdentifier=='+str(atom)+'|| '
        
        return expression
    
    
    
    neighbours = bonddata['neighbours']
    # bondorders = bonddata['bondorders']
    atypes  = bonddata['atypes']
    mtypes     = bonddata['mtypes']
    molecules  = bonddata['molecules']
    
    seeklist    = []
    seekstep    = []
    if dump is not None:
        f = open(dump,'w')
    if shortinfo:
        sf = open(dump[:dump.rfind('.')]+'_shortinfo.txt','w')
    for step,neigh in neighbours.items():
        PRINT = []
        graph = nx.Graph(neigh)
        connected = nx.connected_components(graph)
        for component in connected:
            componentformula = get_molecular_formula(component, atypes, atomsymbols)
            if componentformula == species and component not in seeklist:                
                seeklist.append(component)
                seekstep.append(step)
                # print(step)
                # print(componentformula,component)
                # print(atomids2expression(component))
                if dump is not None:
                    print(step,file=f)
                    print(componentformula,component,file=f)
                    print(atomids2expression(component),file=f)
                if shortinfo:
                    print(step,file=sf)
                    print(componentformula+' =',end=' ',file=sf)
                if source:
                    src = {mtypes[x] for x in component}
                    for j,s in enumerate(src):
                        # print('Source {}:'.format(s))
                        # print(atomids2expression(molecules[s]))
                        # print('-------------------')
                        if dump is not None:
                            print('Source {}:'.format(s),file=f)
                            print(atomids2expression(molecules[s]),file=f)
                            print('-------------------',file=f)
                        if shortinfo:
                            PRINT.append('M-{}'.format(s))
                # print('-'*200)
                if dump is not None: print('-'*200,file=f)
                if shortinfo:
                    print(*PRINT,sep=' + ',file=sf)
                    print(*PRINT,sep=' + ')
    f.close()
    return seeklist,seekstep

@function_runtime
def cumulative_nspecies(neighbours,atypes,atomsymbols,specieslist,**kwargs):
    ###########-get kwargs-###############
    skipts    = kwargs.get('skipts',None)
    ts        = kwargs.get('ts',0.25)
    ######################################
    
    count = [0]*len(specieslist)
    steps = list(neighbours.keys())
    ccount = [[] for i in range(len(specieslist))]
    checklist = [[] for i in range(len(specieslist))]
    for step,neigh in neighbours.items():
        molecules = get_molecules(neigh)
        for molecule in molecules:
            formula = get_molecular_formula(molecule, atypes, atomsymbols)
            for i,species in enumerate(specieslist):
                if species==formula and molecule not in checklist[i]:
                    checklist[i].append(molecule)
                    count[i]+=1
        for i in range(len(ccount)):
            ccount[i].append(count[i])
        
    ps     = step2picosecond(steps, ts)
    
    if skipts is not None:
        index  = ps.index(skipts)
        ps     = ps[index:]
        for i in range(len(ccount)):
            ccount[i] = ccount[i][index:]
    
    ps = np.array(ps)-min(ps)
    print(len(ccount))
    ccount = [np.array(x) for x in ccount]  
    
    return ps,ccount
    
    
    
    
    
    
    
    
    