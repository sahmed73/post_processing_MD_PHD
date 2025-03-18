# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Mar 18 12:21:41 2025

## Purpose:
This script processes bond information from a LAMMPS `bonds.reaxc` file, extracts molecular species, 
and tracks their evolution over time. It allows the user to set a **bond order cutoff** to define connectivity.

"""

import networkx as nx
import pandas as pd

def parsebondfile(bondfilepath, cutoff=0.3, **kwargs):
    """
    Parses the LAMMPS bond output file and extracts molecular connectivity.

    Parameters:
    - bondfilepath (str): Path to the bond output file (`bonds.reaxc`).
    - cutoff (float): Bond order cutoff to determine bonded interactions.
    
    Optional Keyword Arguments:
    - bo (bool): Extract bond orders.
    - mtypes (bool): Extract molecule IDs.
    - charge (bool): Extract atomic charges.
    - abo (bool): Extract atomic bond orders.
    - nlp (bool): Extract number of lone pairs.
    - mols (bool): Extract molecule types.
    - ALL (bool): Extract everything.

    Returns:
    - bonddata (dict): A dictionary containing bond and atomic properties.

    Usage Example:
    bonddata = parsebondfile("bonds.reaxc", cutoff=0.3, bo=True)
    """

    print(f"Bond order cutoff {cutoff} is used")

    # ------- Extract keyword arguments -------------
    bo = kwargs.get('bo', False)  # Extract bond orders
    mtypes = kwargs.get('mtypes', False)  # Extract molecule IDs
    charge = kwargs.get('charge', False)  # Extract atomic charges
    abo = kwargs.get('abo', False)  # Extract atomic bond orders
    nlp = kwargs.get('nlp', False)  # Extract number of lone pairs
    mols = kwargs.get('mols', False)  # Extract molecule types
    ALL = kwargs.get('ALL', False)  # Extract all properties

    def parse():
        """Parses the bond file and extracts atomic connectivity information."""
        bonddata = {}
        neighbours = {}  # Stores connectivity information
        atomtypes = {}  # Stores atom type information
        
        # Conditional data extraction
        if bo or ALL: bondorders = {}
        if mtypes or ALL: molecule_types = {}
        if charge or ALL: charges = {}
        if nlp or ALL: NLP = {}
        if abo or ALL: ABO = {}
        if mols or ALL: molecules = {}

        with open(bondfilepath) as bf:
            prev_natoms = 0
            natoms_flag = False
            warning_flag = False
            first_warning_ignore = True
            atom_counter = 0
            
            for line in bf:
                splitted = line.split()
                
                # Identify new timestep
                if line.find('Timestep') != -1:
                    step = int(splitted[-1])
                    if bo or ALL: bondorders[step] = {}
                    if charge or ALL: charges[step] = {}
                    if nlp or ALL: NLP[step] = {}
                    if abo or ALL: ABO[step] = {}

                    neighbours[step] = {}
                    atomtypes = {}

                # Extract number of particles (atoms)
                if line.find('Number of particles') != -1:
                    current_natoms = int(splitted[-1])

                    if atom_counter != current_natoms:
                        if first_warning_ignore:
                            first_warning_ignore = False
                        else:
                            warning_flag = True

                    if natoms_flag and current_natoms != prev_natoms:
                        print('User warning: Lost atom warning detected!')

                    atom_counter = 0
                    prev_natoms = current_natoms  # Store atom count

                # Process atom information
                if splitted and splitted[0].isnumeric():
                    atom_counter += 1
                    parent = int(splitted[0])  # Atom ID
                    parent_type = int(splitted[1])  # Atom type
                    number_of_children = int(splitted[2])  # Number of bonded atoms
                    children = list(map(int, splitted[3:3 + number_of_children]))  # Bonded atoms

                    # Extract bond orders
                    bo_id1 = 4 + number_of_children
                    bo_idn = bo_id1 + number_of_children
                    bo_values = list(map(float, splitted[bo_id1:bo_idn]))

                    # Apply bond order cutoff
                    updated_children = [child for child, bond_order in zip(children, bo_values) if bond_order >= cutoff]

                    # Store connectivity
                    neighbours[step][parent] = updated_children
                    atomtypes[parent] = parent_type

                    if mtypes or ALL:
                        molecule_types[parent] = int(splitted[3 + number_of_children])
                    if bo or ALL:
                        bondorders[step][parent] = {child: bo_values[i] for i, child in enumerate(children)}
                    if charge or ALL:
                        charges[step][parent] = float(splitted[-1])
                    if nlp or ALL:
                        NLP[step][parent] = float(splitted[-2])
                    if abo or ALL:
                        ABO[step][parent] = float(splitted[-3])

            # Warning for inconsistent atom count
            if warning_flag:
                print('User warning: Repeated atom information detected in some timesteps!')

            # Store extracted data
            bonddata['neighbours'] = neighbours
            bonddata['atypes'] = atomtypes
            if bo or ALL: bonddata['bondorders'] = bondorders
            if mtypes or ALL: bonddata['mtypes'] = molecule_types
            if charge or ALL: bonddata['charge'] = charges
            if nlp or ALL: bonddata['nlp'] = NLP
            if abo or ALL: bonddata['abo'] = ABO
            
            return bonddata

    return parse()

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

def get_species(bondfilepath: str, atomsymbols: list, cutoff: float = 0.3, timestep: float = None):
    """
    Extracts species time series from bond information.

    Parameters:
    - bondfilepath (str): Path to `bonds.reaxc` file from LAMMPS.
    - atomsymbols (list of str): List of element symbols (e.g., ['H', 'C', 'O']).
    - cutoff (float): Bond order cutoff.
    - timestep (float, optional): Time step in femtoseconds (fs), converted to picoseconds (ps).

    Returns:
    - pd.DataFrame: Species count over time.

    Usage:
    df = get_species("bonds.reaxc", ['H', 'C', 'O'], cutoff=0.3, timestep=0.25)
    """
    
    bonddata = parsebondfile(bondfilepath, cutoff=cutoff)
    neighbours = bonddata['neighbours']
    atypes = bonddata['atypes']
    
    count = {}
    for step, neigh in neighbours.items():
        graph = nx.Graph(neigh)
        molecules = nx.connected_components(graph)
        count[step] = {}
        for molecule in molecules:
            species = get_molecular_formula(molecule, atypes, atomsymbols)
            count[step][species] = count[step].get(species, 0) + 1

    df = pd.DataFrame(count).fillna(0).T
    df.index.name = 'Timestep'

    if timestep is not None:
        df['Time'] = df.index * timestep / 1000  # Convert fs to ps
        df.index = df['Time']
        df.drop(['Time'], axis=1, inplace=True)

    return df


#%% User Input Section
dirr = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAO4\25_PAO4_200_O2_Soria\Production\1600\Sim-1"

bondfile_path = dirr + r"\\bonds.reaxc"

# Map atom types to their corresponding element symbols  
# Ensure the order matches the LAMMPS input file  
# Example: pair_coeff * * ${potential_file} H C O  

atomsymbols = ['H', 'C', 'O']


# Get species time series
df = get_species(bondfile_path, atomsymbols, timestep=0.25, cutoff=0.3)

print(df)

#%% Plotting Section
import matplotlib.pyplot as plt

plt.plot(df['H62C30'], label='$H_{62}C_{30}$', color='k')
plt.plot(df['O2'], label='$O_{2}$', color='red')
plt.xlabel('Time (ps)')
plt.ylabel('Number of Species')
plt.legend()
plt.show()
