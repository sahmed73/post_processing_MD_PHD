# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Feb 12 12:42:20 2025
"""

def xyz_to_gaussian_input(xyz_filename, gaussian_filename,
                          method="B3LYP/6-31G*", charge=0, mult=1,
                          nproc=20, mem="90Gb"):
    """
    Converts an XYZ file to a Gaussian input file for DFT optimization.
    
    Parameters:
    xyz_filename (str): Path to the XYZ file.
    gaussian_filename (str): Path to save the Gaussian input file.
    method (str): DFT method and basis set (default: "B3LYP/6-31G*").
    nproc (int): Number of processors (default: 48).
    mem (str): Memory allocation (default: "90Gb").
    
    multiplicity:
        1 (Singlet) → All electrons are paired (e.g., most closed-shell molecules).
        2 (Doublet) → One unpaired electron (e.g., radicals).
        3 (Triplet) → Two unpaired electrons with parallel spins (e.g., oxygen molecule in its ground state).
        4 (Quartet) → Three unpaired electrons, etc.
    """
    
    with open(xyz_filename, 'r') as xyz_file:
        lines = xyz_file.readlines()
    
    # Extract number of atoms and comment line
    num_atoms = int(lines[0].strip())
    molecule_description = lines[1].strip()
    if not molecule_description:
        molecule_description="gaussian optimization"
        print(molecule_description)
    
    # Extract atomic coordinates
    coordinates = lines[2:num_atoms+2]
    
    # Write Gaussian input file
    with open(gaussian_filename, 'w') as gaussian_file:
        gaussian_file.write(f"%NProcShared={nproc}\n")
        gaussian_file.write(f"%mem={mem}\n")
        gaussian_file.write(f"#n {method} Opt freq\n\n")
        gaussian_file.write(f"{molecule_description}\n\n")
        gaussian_file.write(f"{charge} {mult}\n")  # Default charge and multiplicity (neutral singlet)
        
        for line in coordinates:
            gaussian_file.write(line)
        
        gaussian_file.write("\n\n")
    
    print(f"Gaussian input file saved as {gaussian_filename}")

def creat_submit_file(submit_filename, job_name, partition):
    content=f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --partition {partition}
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=20
#SBATCH --export=ALL
#SBATCH --mail-user=shihab2196@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --mem=90G

export SLURM_CLUSTERS=merced

inname=gaussian_opt.in
jobname=$(basename -s .in $inname)
outname=$jobname.out

module load gaussian/g16-b01
#module load gaussian/g09-d01

g16 $inname $outname'''
    
    with open(submit_filename, 'w') as sfile:
        sfile.write(content)



# Example usage
which =     1
dirr = rf"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\GAUSSIAN\OPT\JOB-0002_A1-A5_OPT+Freq\A0004\post{which}"
xyzfile = dirr + rf"\A0004_un-optimized_post{which}.xyz"
gaussian_input_filename = dirr+r"\gaussian_opt.in"
xyz_to_gaussian_input(xyzfile, gaussian_input_filename, charge=0, mult=2)
creat_submit_file(dirr+"\\submit.sh", job_name=f"A4_post{which}_opt",
                  partition="pi.amartini")
