# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Tue Nov  5 12:43:51 2024
"""

import re

def extract_lammps_variables(input_file_path):
    """
    Extracts variables defined in a LAMMPS input file.

    This function reads a LAMMPS input file and identifies variables defined using the 
    `variable` command. It supports all standard LAMMPS variable styles, including:
    `delete`, `atomfile`, `file`, `format`, `getenv`, `index`, `internal`, `loop`,
    `python`, `string`, `timer`, `uloop`, `universe`, `world`, `equal`, `vector`, and `atom`.

    For each variable, the function captures:
      - The variable name
      - The variable style (type of variable)
      - The variable value or expression

    Inline comments, which follow a `#` symbol, are ignored. The function also skips
    any lines that are fully commented out.

    Parameters:
        input_file_path (str): The path to the LAMMPS input file to be read.

    Returns:
        dict: A dictionary containing variable definitions where each key is a variable name,
              and the value is another dictionary with the following keys:
              - 'style' (str): The variable style (e.g., `equal`, `string`, etc.).
              - 'value' (str): The assigned value or expression for the variable.

    Example:
        Given an input file containing:
            variable temp equal 300 # temperature in kelvin
            variable density string "1.0 g/cm^3" # density of the system

        The function returns:
            {
                'temp': {'style': 'equal', 'value': '300'},
                'density': {'style': 'string', 'value': '1.0 g/cm^3'}
            }
    """
    
    variables = {}

    # Regular expression to match the LAMMPS `variable` command with all styles
    styles = "(delete|atomfile|file|format|getenv|index|internal|loop|python|string|timer|uloop|universe|world|equal|vector|atom)"
    variable_pattern = re.compile(rf"^\s*variable\s+(\w+)\s+{styles}\s+(.*)", re.IGNORECASE)

    # Open and read the LAMMPS input file
    with open(input_file_path, 'r') as file:
        for line in file:
            # Remove any inline comment by splitting at the first occurrence of #
            line = line.split('#', 1)[0].strip()
            
            # Skip empty lines (lines that were only comments)
            if not line:
                continue
            
            # Look for lines starting with 'variable' followed by a valid style
            match = variable_pattern.match(line)
            if match:
                var_name = match.group(1)  # Variable name
                var_style = match.group(2)  # Variable style (delete, file, format, etc.)
                var_value = match.group(3).strip()  # Variable value or expression
                
                # Store in dictionary with variable name as key
                variables[var_name] = {'style': var_style, 'value': var_value}
    
    return variables


input_file_path = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\ARLOI\ARLOI-V1\11_ARLOI_Slow_Ramp_1KpPS_300-500K\A0004\Equilibration\NVT-2\input.in"  # Replace with your LAMMPS input file path
variables = extract_lammps_variables(input_file_path)

# Print the extracted variables
for var_name, var_info in variables.items():
    print(f"Variable: {var_name}, Style: {var_info['style']}, Value: {var_info['value']}")
