# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 00:50:24 2023

@author: arup2
"""

import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(nrows=50, ncols=1, figsize=(12, 12))
n  = 3
nn = 100
for i in range(len(ax)):
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    ax[i].set_xlim(-0.1,n+0.1)
    [ax[i].spines[x].set_visible(False) for x in ['top','right','left','bottom']]
    ax[i].plot([3/(50-i+1)]*nn,np.linspace(0,n,nn),linewidth=5,color='maroon')
    ax[i].plot(np.linspace(0,n,nn),[n/2]*nn,color='k')

[ax[i].spines[x].set_visible(True) for x in ['bottom']]
ax[i].set_xticks(range(0,n+1))