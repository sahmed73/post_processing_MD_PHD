# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 31 13:52:47 2024
"""

import os
import re
import sys

def calculate_total_cpu_hours_from_log(file_path):
    total_cpu_hours = 0.0
    pattern = r"Loop time of ([\d.]+) on (\d+) procs"

    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                loop_time = float(match.group(1))
                num_procs = int(match.group(2))
                total_cpu_hours += (loop_time / 3600) * num_procs

    return total_cpu_hours

def calculate_total_cpu_hours_in_directory(root_dir):
    total_cpu_hours = 0.0
    log_file_count = 0

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename == "log.lammps":
                log_path = os.path.join(dirpath, filename)
                total_cpu_hours += calculate_total_cpu_hours_from_log(log_path)
                log_file_count += 1

    return total_cpu_hours, log_file_count

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_cpu_hours.py <root_directory>")
        sys.exit(1)

    root_directory = sys.argv[1]
    total_cpu_hours, log_file_count = calculate_total_cpu_hours_in_directory(root_directory)
    print("Total CPU hours:", total_cpu_hours)
    print("Number of log files encountered:", log_file_count)

