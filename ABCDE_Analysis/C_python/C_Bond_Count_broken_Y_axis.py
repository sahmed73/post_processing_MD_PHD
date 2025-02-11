# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 03:20:08 2023

@author: arup2
"""
import magnolia.bondfile_parser as bfp
import matplotlib.pyplot as plt
import magnolia.needless_essential as ne

directory = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\C\C_300_O2\Production\LessThan_000K_1400'
filename  = '\\bonds.reaxc'
bondfile  = directory+filename

atomsymbols = ['H','C','O']
bondinfo = bfp.parsebondfile(bondfile,bo=True)
#%%--
neighbours = bondinfo['neighbours']
atypes     = bondinfo['atypes']
bondorders = bondinfo['bondorders']
steps      = list(neighbours.keys())

####-Bond Extracting-##############
firststepneigh = neighbours[steps[0]]
# 1 H
# 2 C
# 3 O
bonds = {}

CH3_O = [] # to find bond-5
for parent,children in firststepneigh.items():
    # O-H bonds    
    if atypes[parent] == 3:                
        for child in children:
            if atypes[child]==1:
                key = 'O-H'
                if key not in bonds.keys():
                    bonds[key] = [(parent,child)]
                else:
                    bonds[key].append((parent,child))
    
    # C-C_link bonds
    if atypes[parent] == 2: # Carbon
        H_count = [atypes[x] for x in children].count(1) # Hydrogen
        C_count = [atypes[x] for x in children].count(2) # Carbon
        if H_count == 2 and C_count == 2:
            for child in children:
                if atypes[child]==2:
                    key = 'C-C_link'
                    if key not in bonds.keys():
                        bonds[key] = [(parent,child)]
                    else:
                        bonds[key].append((parent,child))
    
    # Tertiary Butyl link bonds
    if atypes[parent] == 2: #Carbon
        key = 'tertiary'
        if [atypes[x] for x in children].count(2)==4:
            for child in children:
                if [atypes[x] for x in firststepneigh[child]].count(2)==3:
                    if key not in bonds.keys():
                        bonds[key] = [(parent,child)]
                    else:
                        bonds[key].append((parent,child))
                        
                    ## -CH3 Bonds
                    for carbon in [parent,child]:
                        ff = firststepneigh[carbon]
                        if [atypes[x] for x in ff].count(2)==4:
                            ttchildren = firststepneigh[carbon]
                            for ttchild in ttchildren:
                                keyy = '-CH3'
                                fff = firststepneigh[ttchild]
                                if [atypes[x] for x in fff].count(1)==3:
                                    CH3_O.append(ttchild)
                                    if keyy not in bonds.keys():
                                        bonds[keyy] = [(carbon,ttchild)]
                                    else:
                                        bonds[keyy].append((carbon,ttchild))
### finding extra bonds #######
key  = 'H2C=O'
keyy = 'O-H formation'
bonds[key]  = []
bonds[keyy] = []
for step, neigh in neighbours.items():
    for parent, children in neigh.items():
        
        # finding H2C=O
        Noxy = [atypes[x] for x in children].count(3)
        Nhyd = [atypes[x] for x in children].count(1)
        if parent in CH3_O and Noxy==1 and Nhyd==2:
            for child in children:
                bond = (parent,child)
                bo = bondorders[step][parent][child]
                if child>3750 and bo>1.4 and bond not in bonds[key]:
                    print(bo)
                    bonds[key].append(bond)
        
        # finding O-H fromation
        # if parent>3750:
        #     for child in children:
        #         if len(children)==2 and [atypes[x] for x in children].count(1)==2:
        #             bond = (parent,child)
        #             bonds[keyy].append(bond)

bonds[keyy] = [16*x for x in bonds[key]]
               
bondname = {'O-H':'Bond-1 dissociation','C-C_link':'Bond-2 dissociation',
            'tertiary':'Bond-3 dissociation', '-CH3':'Bond-4 dissociation',
            'H2C=O':'(H$_2$C)=O bond formation',
            'O-H formation':'O-H bond formation'}
for key in bondname:
    print(key,'--',len(bondname[key]))
    

# adding missing O-H bond
missing = (256,300)
bonds['O-H'].append(missing)
print('added missing O-H bond',missing)

# delete extra C-C_link bonds
extra = [(3301,3302),(3301,3306),(3302,3301),(3306,3301)]
bond = bonds['C-C_link']
bond = [x for x in bond if x not in extra]
bonds['C-C_link'] = bond

skipts = 70
savedir = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\BO_Evolution'
molecule = directory[directory.find('ABCDE')+8]
print(molecule)
title  = directory[directory.find('ABCDE'):]+'\nMolecule {}: {} bonds'.format(molecule,key)
#%%-
bond_type = list(bondname.keys())
marker = ['s','^','x','o','>','*']    
color    = ['b','c','m','r','g','grey']

cutoff = 1
tol    = 0.3

plot_keys = bond_type
print(plot_keys)
print(bonds.keys())

fig, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.15)

for i,key in enumerate(plot_keys):
    if key in plot_keys:
        bondlist = bonds[key]
        markevery = 150
        
        if key=='O-H formation':
            bondlist = bonds['H2C=O']
            ps, nbonds = bfp.get_nbondsVStime(bondorders, atypes,
                                          bondlist,skipts=skipts,cutoff=cutoff,
                                          tol=tol)
            nbonds = nbonds*9
        else:
            ps, nbonds = bfp.get_nbondsVStime(bondorders, atypes,
                                          bondlist,skipts=skipts,cutoff=cutoff,
                                          tol=tol)
        
        ax1.plot(ps,nbonds,label=bondname[key],marker=marker[i],
                 markevery=markevery,markersize=4,linewidth=0.8,color=color[i])
        ax2.plot(ps,nbonds,label=bondname[key],marker=marker[i],
                 markevery=markevery,markersize=4,linewidth=0.8,color=color[i])

###-plot-####

ax1.set_ylim(0, 201)  # outliers only
ax2.set_yticks(range(450,650+1,50))
ax2.set_ylim(450, 650)  # most of the data

# hide the spines between ax and ax2
ax2.spines.bottom.set_visible(False)
ax1.spines.top.set_visible(False)
ax2.xaxis.tick_top()
ax2.tick_params(labeltop=False)  # don't put tick labels at the top
ax1.xaxis.tick_bottom()

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax2.plot([0, 1], [0, 0], transform=ax2.transAxes, **kwargs)
ax1.plot([0, 1], [1, 1], transform=ax1.transAxes, **kwargs)

fig.text(0.03, 0.5, 'Number of bonds', va='center', rotation='vertical',
         fontsize=12)

d = directory[directory.find('LAMMPS')+7:]
title = d[d.find('\\')+1:]+'\n'+'Molecule C\n'
plt.xlim(0,1000)
# plt.title(title)
#plt.title('Number of bond vs time {}'.format(bond_type))
ax2.legend(bbox_to_anchor=(1.51,1))
plt.xlabel('Time (ps)',fontsize=12)
# plt.xlim(right=1150)
plt.savefig('..\\..\\python_outputs\\bondplot\\bondplot-'+ne.randstr(), dpi=300,bbox_inches='tight')
