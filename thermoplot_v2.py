# -*- coding: utf-8 -*-
"""
Created on Sun May 28 01:45:11 2023

@author: arup2
"""

import magnolia.log_parser_FHB as lfp
import matplotlib.pyplot as plt
import numpy as np
import magnolia.needless_essential as ne
plt.style.use('default')
plt.rcParams['font.size'] = 18

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\PAOtcr\20PAOtcr_15A\Equilibrate\200NVT'

filename = '\\log.lammps'
logfile = directory+filename
timestep = 0.25
thermo = lfp.thermo_dict_v2(logfile,serial=2)

plot_property = 'PotEng'
steps         = np.array(thermo['Step'])
ps            = steps*timestep/1000
    

thermo_property = np.array(thermo[plot_property])
if plot_property in ['PotEng','TotEng','KinEng']:
    thermo_property = thermo_property - min(thermo_property)
 
start = 0
end = None

ps = ps-ps[0]
plt.plot(ps[start:end],thermo_property[start:end],color='tab:red')

plt.minorticks_on()
plt.xlabel('Time (ps)') 
plt.ylabel(ne.getlab(plot_property))
# plt.title(directory[directory.find('Base'):]+'\n\n', fontsize=12)
plt.savefig(r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\thermoplot\thermoplot-'+ne.randstr()+'.png', dpi=300,bbox_inches='tight',transparent=True)
print(ps[-1])
plt.show()