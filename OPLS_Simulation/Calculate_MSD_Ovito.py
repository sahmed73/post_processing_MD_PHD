# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Apr  8 00:23:32 2024
"""

from ovito.io import import_file, export_file
from ovito.modifiers import CalculateDisplacementsModifier, ExpressionSelectionModifier, DeleteSelectedModifier
import numpy as np
import matplotlib.pyplot as plt
import magnolia.plot_template as mplt

def calculate_msd(frame, data):
    displacement_magnitudes = data.particles['Displacement Magnitude']
    msd = np.sum(displacement_magnitudes ** 2) / len(displacement_magnitudes)
    data.attributes["MSD"] = msd 

antioxidants = ['A', 'B', 'BHT']
temperatures = [500]
nevery = 2
markevery = 50

colors = ['tab:blue', 'tab:red', 'tab:green']
markers  = ['o', 's', '<']

for i, ao in enumerate(antioxidants):
    for j, temp in enumerate(temperatures):
        print(ao,temp)
        pipeline = import_file(r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\OPLSAA\PAO_Antioxidant\PAO+{}\10PAO+3{}\Single\{}K\Sim-1\production.nvt.dump".format(ao,ao,temp))
        
        # Calculate per-particle displacements with respect to initial simulation frame:
        pipeline.modifiers.append(CalculateDisplacementsModifier())
        
        # select all atoms that are NOT antioxidants
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='MoleculeIdentifier<11'))
        
        # Delete all selected atoms, leaving only those belonging to the specified molecule.
        pipeline.modifiers.append(DeleteSelectedModifier())        
        pipeline.modifiers.append(calculate_msd)
        
        msd_data = []
        steps    = []
        nframes = pipeline.source.num_frames
        for frame in range(nframes):
            data = pipeline.compute(frame)
            msd_data.append(data.attributes['MSD'])
            steps.append(data.attributes['Timestep'])
            
        msd_data = np.array(msd_data)
        steps = np.array(steps)*0.25/1e6
        plt.plot(steps[::nevery],msd_data[::nevery],c=colors[i],
                 label=f'{ao}')
        
        
plt.xlabel('Time (ns)')
plt.ylabel('MSD')
plt.title('{}K\n'.format(*temperatures))
plt.legend(loc='upper left')
plt.ylim(top=60)
mplt.saveplot(name='MSD_A_B_BHT_{}'.format(*temperatures))
