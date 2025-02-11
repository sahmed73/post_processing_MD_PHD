# -*- coding: utf-8 -*-
"""
Created on Fri May 19 04:40:03 2023

@author: arup2
"""
import matplotlib.pyplot as plt

onsets  = [1351,1104,1513]
error   = [135.5,47.9,62.4]
lignins = ['B','D','E']
xvals   = [1,2,3]

plt.style.use('classic')
light_maroon = 'maroon'
plt.bar(lignins,onsets,yerr=error,color=light_maroon,edgecolor='black',capsize=7,width = 0.45,align='center',error_kw=dict(lw=2, capthick=2))
#plt.xticks(xvals, lignins)
plt.xlim([- 0.5, 2 + 0.5])

plt.ylabel('Onset temperature (K)')
#plt.xlabel('Molecule')
plt.savefig('python_outputs\\bargraph_onset_lignin.png',dpi=400)