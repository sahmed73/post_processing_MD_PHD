# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Aug 26 13:22:00 2024
"""
import numpy as np
import matplotlib.pyplot as plt

h = np.array([0.25,0.28475])*4*1000+300
src = np.array([7.928,11.82])
x = ['BHT', 'A004']

# Define colors and shapes for each point
colors = ['red', 'blue', 'green', 'purple', 'k']
shapes = ['o', 's', '^', 'D', '>']  # 'o' for circle, 's' for square, '^' for triangle, 'D' for diamond

# Create scatter plot
plt.figure(figsize=(8, 6))

for i in range(h.size):
    plt.scatter(x[i], src[i], color=colors[i], marker=shapes[i], s=200, label=f'Point {i+1}')

plt.xlabel('')
# plt.ylim(top=1500+10)
plt.ylabel('SRC ($ns^{-1}$)')
# plt.legend()
plt.show()