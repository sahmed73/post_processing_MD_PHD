# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 15 10:40:22 2024
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Creating the data

RSI=np.array([5.20693664222082, 1.7369663301369145, 4.492157287325462, 4.555204057497001, 7.742859409407634])
N_OH=np.array([1,2,1,1,2])
data = {
    'Antioxidant': ['A0001', 'A0002', 'A0003', 'A0004', 'A0005'],
    'RSI per Ph-OH': RSI/N_OH
}

# Creating the DataFrame
df_bar = pd.DataFrame(data)

# Plotting a simple bar plot
fig, ax  = plt.subplots(dpi=350)
sns.barplot(data=df_bar, x='Antioxidant', y='RSI per Ph-OH', palette="tab10", ax=ax)
