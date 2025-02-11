# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 14 22:43:58 2024
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define Antioxidants (AOs)
AOs = [f'A{i+1:04}' for i in range(5)]
print(AOs)

# Initialize an empty DataFrame
df = pd.DataFrame()

# Loop through each AO and aggregate data
for AO in AOs:
    dirr = rf'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\{AO}\Production\300-500K_TRR=1Kpps'
    infile = dirr + r'\scavenged_H_counts.csv'
    dff = pd.read_csv(infile, index_col=0)
    dff['Antioxidant'] = AO
    df = pd.concat([df, dff])

# Convert to long format for Seaborn

df_long = df.reset_index().melt(id_vars=['index', 'Antioxidant'], var_name='Source', value_name='Counts')

# Define your custom order for the x-axis labels
custom_order=['Ph—OH', 'Self', 't-Bu', 'Others', 'non-Ph—OH', 'Ar']  # Replace with your specific order

# Set the custom order for the 'index' column
df_long['index'] = pd.Categorical(df_long['index'], categories=custom_order, ordered=True)

# Plotting
plt.rcParams['font.size'] = 15
fig, ax = plt.subplots(dpi=350, figsize=[6.4, 4.8])

# Let Seaborn calculate the mean and 95% CI
sns.barplot(data=df_long, x='index', y='Counts', hue='Antioxidant',
            errorbar=('ci', 95), capsize=0.2, ax=ax)

# Customize plot
ax.set_xlabel("Source of Scavenger H")
ax.set_ylabel("Average Counts")
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()


