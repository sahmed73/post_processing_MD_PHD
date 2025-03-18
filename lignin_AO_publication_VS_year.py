# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Mar  7 14:10:51 2025
"""

import matplotlib.pyplot as plt
import numpy as np

# Data from Figure 22
years = np.array(range(2008, 2023))
articles = np.array([24, 30, 38, 45, 54, 65, 74, 81, 110, 139, 153, 202, 239, 341, 302])

# Create a cleaner plot
fig, ax = plt.subplots(figsize=(8, 5), dpi=350)

ax.plot(years, articles, marker='o', linestyle='-', color='grey', markerfacecolor='black', linewidth=2)

# Improve x-axis labels
ax.set_xticks(years)
ax.set_xticklabels(years, rotation=45)
ax.tick_params(axis='both', labelsize=14)

# Improve y-axis labels
ax.set_yticks(np.arange(0, 401, 50))

# Labels and title
ax.set_xlabel("Year", fontsize=20, labelpad=10)
ax.set_ylabel("Number of Articles", fontsize=20, labelpad=10)
ax.set_title("Trend of Publications on Lignin Antioxidants", fontsize=16, pad=15)

# Grid and aesthetics
ax.grid(True, linestyle="--", alpha=0.4)
ax.set_facecolor("lightblue")  # Light background for clarity
fig.patch.set_facecolor("white")

# Show the improved plot
plt.show()
