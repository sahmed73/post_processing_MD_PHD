# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Oct 23 05:17:08 2023
"""

import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import sys
import magnolia.needless_essential as ne
import pandas as pd
import numpy as np
import seaborn as sns
import random
import networkx as nx


directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\Squalane\25_Squalane_200_O2_Soria\Production\1600\Sim-1'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

asyms = ['H','C','O']
atominfo = bfp.get_neighbours(bondfile,bo=True)
#%%
neighbours, atypes, bondorders = atominfo

firststepneigh = neighbours[list(neighbours.keys())[0]]
graph = nx.Graph(firststepneigh)

sym_edges = []

for i in range(0,25):
    nodes = range(1+92*i,92+92*i+1)
    edges = graph.subgraph(nodes).edges
    sym_edge = [(atypes[x],atypes[y]) for x,y in edges]
    sym_edges.append(sym_edge)

array = np.array(sym_edges)
transposed_arr = array.transpose(1, 0, 2)
lst = transposed_arr.tolist()

for a in lst:
    first = a[0]
    # print(a.count(first))
    print(a)
    print('--------------')