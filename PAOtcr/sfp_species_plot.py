# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Apr 18 01:54:20 2024
"""

import magnolia.speciesfile_parser as sfp
import magnolia.bondfile_parser as bfp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt
mplt.custom_plot_features(fontsize_incr=2)

fig, axs = plt.subplots(1,3, sharey=True, figsize=[5,3], dpi=300)
fig.subplots_adjust(wspace=0.5)
for k, temp in enumerate([600,700,800]):
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\3PAOtcr_3A_Using_PACKMOL\Production\StartedFromEnergyEquilibration\{}K'.format(temp)
    filename = '\\species.out'
    speciesfile = dirr+filename
    species = sfp.get_species_count(speciesfile, timestep=0.25)
    colors = ['b','r']
    for i, spec in enumerate(species):
        if spec=='Time': continue
        if spec not in ['H61C30O','H30C19O3']: continue
        label = bfp.make_molecular_formula_latex(spec,sort=True)
        axs[k].scatter(species['Time']/1000,species[spec], s=1,
                    edgecolors=colors[i], alpha=0.1)
        axs[k].set_xlim(0,1)
        axs[k].grid(False)
        axs[k].set_ylim(-0.1,3.1)
        axs[k].set_xticks([0,0.5,1.0])
        axs[k].set_yticks([0,1,2,3])
        axs[k].tick_params(axis='y', which='major', length=0)
        if k==len(axs)-1:
            axs[k].plot([],[],label=label,c=colors[i],linewidth=3)

axs[k].legend(loc=[1.05,0.7])

fig.text(0.5, -0.06, 'Time (ns)', ha='center')
axs[0].set_ylabel('Number of molecule')
mplt.saveplot('species','species')