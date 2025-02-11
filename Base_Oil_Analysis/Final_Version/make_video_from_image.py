# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Jan 26 01:17:27 2024
"""

import os
import imageio

# Path to the GIF file
folder_path = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs\pathway'

# Output GIF file name (e.g., output.gif)
output_video = folder_path+'\\output.mp4'

# Create a list of PNG file paths sorted by filename
png_files = sorted([os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.png')],
                   key=lambda x: x[-10:])

# Create a list to store the frames
frames = []

# Read each PNG file and append it to the frames list
for png_file in png_files:
    img = imageio.imread(png_file)
    frames.append(img)

# Create the video using imageio
imageio.mimsave(output_video, frames, fps=1.3)  # Adjust the fps (frames per second) as needed

print(f"Video '{output_video}' created successfully.")



