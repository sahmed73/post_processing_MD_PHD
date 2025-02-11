# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri May 31 09:00:39 2024
"""

import pandas as pd
import sys
import numpy as np
from periodictable import elements

def pdb2lmpdata(pdbfile, atomtype_order=None, box_size=None, box_incr=None,
                adjust_coord=True):
    # atom_types is a dictionary key=symbol, value=atom_id
    # box_size=[0., 48., 0., 48., 0., 48.] # xlo, xhi, ylo, yhi, zlo, zhi
    # if box_size is None then it will assign the min, max coord as dimension
    # box_incr=[0., 2. , 0., 2.,  0. , 2.] # example
    # adjust_coord True: atom coord trans so that min atom pos change to 0
    
    '''
    COLUMNS        DATA TYPE       FIELD         DEFINITION (short)
    ----------------------------------------------------------------------------
      1 - 6        Record name     "ATOM  "
      7 -11        Integer         Atom serial number (atom_id)
    13 - 16        Atom            Atom name (atom_name)
    17             Character       Alternate location indicator
    18 - 20        Residue name    Residue name (res_name)
    22             Character       Chain identifier
    23 - 26        Integer         Residue sequence number (res_id)
    27             AChar           Code for insertion of residues
    31 - 38        Real(8.3)       Orthogonal coordinates for X (X)
    39 - 46        Real(8.3)       Orthogonal coordinates for Y (Y)
    47 - 54        Real(8.3)       Orthogonal coordinates for Z (Z)
    55 - 60        Real(6.2)       Occupancy
    61 - 66        Real(6.2)       Temperature factor
    77 - 78        LString(2)      Element symbol (sym)
    79 - 80        LString(2)      Charge on the atom (charge)
    '''
    
    start = None
    with open(pdbfile, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                end = i
                if start is None:
                    start = i
    
    colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22),
                (22, 26), (26, 27), (30, 38), (38, 46), (46, 54), (54, 60),
                (60, 66), (76, 78), (78, 80)]
    
    names = ['Record name', 'atom_id', 'atom_name',
             'Alternate location indicator',
             'res_name', 'Chain identifier', 'res_id',
             'Insertion code', 'X', 'Y', 'Z', 'Occupancy',
             'Temperature factor', 'sym', 'charge']

    # Read the PDB file
    df_pdb = pd.read_fwf(pdbfile, colspecs=colspecs, header=None, names=names,
                     skiprows=start, nrows=end-start+1)
    
    lmpdata_header = ['atom_id', 'res_id', 'sym', 'charge', 'X', 'Y', 'Z']
    df_lmpdata = df_pdb[lmpdata_header].fillna(0)
    
    # Cahrge handling: convert 1- charges to -1
    def convert_charge(charge):
        if isinstance(charge, str) and charge[-1] in ['-', '+']:
            return charge[::-1]
        return charge
    df_lmpdata['charge'] = df_lmpdata['charge'].apply(convert_charge)
    
    # atom type mapping
    if atomtype_order is None: # alphabatcal order
        atomtype_order = sorted(df_lmpdata['sym'].unique())
    
    atom_type_map = dict(zip(atomtype_order,range(1,1+len(atomtype_order))))
    
    
    df_lmpdata['sym'] = df_lmpdata['sym'].map(atom_type_map)
    df_lmpdata.rename(columns={'sym': 'atom_type'}, inplace=True)
    
    if box_size is None:
        xlo, ylo, zlo = np.round(df_lmpdata[['X','Y','Z']].min())
        xhi, yhi, zhi = np.round(df_lmpdata[['X','Y','Z']].max())
        if box_incr is not None:
            temp_box = np.array([xlo, xhi, ylo, yhi, zlo, zhi])
            xlo, xhi, ylo, yhi, zlo, zhi = temp_box+np.array(box_incr)
    else:
        if box_incr is None:  
            xlo, xhi, ylo, yhi, zlo, zhi = box_size
        else:
            xlo, xhi, ylo, yhi, zlo, zhi = [x+y for x,y in zip(box_size,box_incr)]
            
    # translate the atoms to the middle
    if box_incr is not None:
        xmin, xmax, ymin, ymax, zmin, zmax = box_incr
        trans = np.array([(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2])
        df_lmpdata[['X', 'Y', 'Z']] = df_lmpdata[['X', 'Y', 'Z']] + trans
        
    # adjust coordinate: make min coordinate 0 using translation
    if adjust_coord:
        df_lmpdata['X']-=xlo
        df_lmpdata['Y']-=ylo
        df_lmpdata['Z']-=zlo
        
        temp_box = np.array([xlo, xhi, ylo, yhi, zlo, zhi])
        trans    = np.array([xlo,xlo,ylo,ylo,zlo,zlo])
        # update xlo, xhi, ylo, yhi, zlo, zhi
        xlo, xhi, ylo, yhi, zlo, zhi = temp_box-trans
    
    n_atom_types = len(atom_type_map)
    natoms, _    = df_lmpdata.shape
    
    header = f'''# LAMMPS data file converted PDB file # Created by SA
    
{natoms}  atoms
{n_atom_types}  atom types

{xlo}  {xhi}   xlo xhi
{ylo}  {yhi}   ylo yhi
{zlo}  {zhi}   zlo zhi

Masses

'''
    masses = ''
    for symbol, atom_type in atom_type_map.items():
        element = elements.symbol(symbol)
        mass = f'{atom_type} {element.mass} # {element}\n'
        masses+=mass
    
    print(header+masses)
    output = pdbfile[:pdbfile.rfind('.')]+'.data' 
    string = df_lmpdata.to_string(header=None, index=None)
    with open(output, 'w') as out:
        out.write(header)
        out.write(masses+'\n')
        out.write('Atoms\n\n')
        out.write(string)
    return df_pdb

pdbfile=sys.argv[1]

atomtype_order = ['H','C','O']
box_size       = None #[0., 48., 0., 48., 0., 48.] # xlo, xhi, ylo, yhi, zlo, zhi
box_incr       = [0.,  2., 0.,  2., 0.,  2.]

df = pdb2lmpdata(pdbfile, atomtype_order=atomtype_order,
                     box_size=box_size, box_incr=box_incr)