# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Oct 30 23:17:20 2024
"""

import re

def calculate_total_cpu_hours_from_log(file_path):
    """
    Calculate total CPU hours used across multiple LAMMPS runs by reading a log file.

    Parameters:
    file_path (str): Path to the LAMMPS output or log file.

    Returns:
    float: Total CPU hours across all runs.
    """
    total_cpu_hours = 0.0
    
    # Define a regular expression pattern to capture loop time and processor count
    pattern = r"Loop time of ([\d.]+) on (\d+) procs"

    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                loop_time = float(match.group(1))  # Extract loop time in seconds
                num_procs = int(match.group(2))   # Extract number of processors
                # Calculate CPU hours for this run and add to the total
                total_cpu_hours += (loop_time / 3600) * num_procs

    return total_cpu_hours

# Example usage
file_path = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\10_ARLOI_NPT-high_NPT-cooling_NPT-low_SORIA\A0001\Production\Sim-1\log.lammps"
total_cpu_hours = calculate_total_cpu_hours_from_log(file_path)
print("Total CPU hours:", total_cpu_hours)
