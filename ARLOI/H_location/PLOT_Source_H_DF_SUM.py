# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 14 12:19:37 2024
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

AOs = [f'A{i+1:04}' for i in range(5)]
print(AOs)
df = pd.DataFrame()
for AO in AOs:
    dirr=rf'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\{AO}\Production\300-500K_TRR=1Kpps'
    infile=dirr+r'\scavenged_H_counts.csv'
    dff=pd.read_csv(infile,index_col=0)
    series=dff.sum(axis='columns')
    df[AO]=series
    print(AO)
    print(100*dff.sum(axis=1)/dff.sum().sum())

df_sort=df.sort_values(by='A0001', ascending=False)
df_long = df_sort.reset_index().melt(id_vars='index', var_name='Antioxidant', value_name='Counts')

plt.rcParams['font.size']=15
fig, ax = plt.subplots(dpi=350, figsize=[6.0, 4.0])
sns.barplot(data=df_long, x='index', y='Counts', hue='Antioxidant', ax=ax)
ax.set_xlabel("Source of Scavenger H")
plt.xticks(rotation=45)
