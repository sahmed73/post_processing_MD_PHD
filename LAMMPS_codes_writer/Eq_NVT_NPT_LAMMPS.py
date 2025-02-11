# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Apr  1 01:00:07 2024
"""
import tkinter as tk
from tkinter import simpledialog
import os

class CustomMultiInputDialog(simpledialog.Dialog):
    def __init__(self, master, inputs):
        self.inputs = inputs
        super().__init__(master, title="Input Parameters")

    def body(self, master):
        self.entries = {}
        for i, (label, prefill) in enumerate(self.inputs.items()):
            tk.Label(master, text=f"{label}:").grid(row=i, column=0, sticky='w')
            entry = tk.Entry(master, width=50)
            entry.grid(row=i, column=1)
            entry.insert(0, "" if prefill is None else str(prefill))
            self.entries[label] = entry

        return list(self.entries.values())[0]  # Set focus to the first entry widget

    def apply(self):
        # Store all entry values in a dictionary, converting strings back to original types where necessary
        self.result = {label: self._convert_type(entry.get()) for label, entry in self.entries.items()}

    def _convert_type(self, value):
        # Convert value back to original type if possible, for 'None', int, float, and str
        if value == "":
            return None
        try:
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return value

inputs = dict(
    T0               = 300,
    T                = 300,
    P0               = 1,
    NVT_runtime      = 200,   # ps
    NPT_runtime      = 1000,  # ps
    tstep            = 0.25,  # fs
    datafile         = "../../DataFile/20_PAO_Radical_20_BHT.data",
    restart_file     = None,
    res_filename     = "repeated_eq_restart",
    res_freq         = 10,
    sim_name         = "equilibrated",
    sim_number       = 1,
    potential_file   = "../../../../Potentials/CHONSSi_Soria.ff",
    thermostep       = 1000,
    dumpstep         = 1000,
    
    ######## submit file info #######
    job_name         = "Eq-20rPAO-20BHT-Soria",
    partition        = "long",
    partition_merced = "compute",
    nodes            = 1,
    nodes_merced     = 1,
    ntask            = 36,
    ntask_merced     = 32,
    time             = "2-23:59:00",
    time_merced      = "119:59:00"
)

# Create the application window
root = tk.Tk()
root.withdraw()  # Hide the main window

# Create and show the dialog, and wait for the user to close it
dialog = CustomMultiInputDialog(root, inputs)
user_inputs = dialog.result  # After the dialog is closed, get the result

# Use the input values
# print("User entered values:")
# for key, value in user_inputs.items():
#     print(f"{key}: {value}")



initialize_variables = '''#Equilibration Using NVT and NPT
#Initialize variables--------------------------------------------------------
variable			T0 equal {T0}
variable			T equal {T}
variable			P0 equal {P0}
variable			runtime equal {NVT_runtime}
variable			tstep equal {tstep}
variable			datafile string {datafile}
variable			restart_file string {restart_file}
variable			res_filename string {res_filename}
variable			res_freq equal {res_freq}
variable			sim_name string {sim_name}
variable			sim_number equal {sim_number}
variable			potential_file string {potential_file}
variable			thermostep equal {thermostep}
variable			dumpstep equal {dumpstep}
'''.format(**user_inputs)

body_nvt = '''
#Define Variable ------------------------------------------------------------
variable			runstep equal 1000*${runtime}/${tstep}
variable			resstep equal ${runstep}/${res_freq}
variable 		seed equal (289567+87659*${sim_number})

#Set up simulation-----------------------------------------------------------
units				real
boundary			p p p
atom_style			full/kk
neighbor			2.0 bin
neigh_modify		delay 0 every 10 check no
timestep			${tstep}
read_data			${datafile}
reset_timestep		0
restart				${resstep} ${res_filename}

#Define potential information------------------------------------------------
pair_style         reaxff/kk NULL checkqeq yes
pair_coeff         * * ${potential_file} H C O

#Include mandatory ReaxFF charge equilibration (QEq) fix--------------------
fix		   		   0 all qeq/reaxff/kk 1 0.0 10.0 1.0e-6 reaxff

#Minimize the unit cell--------------------------------------------------------
min_style          cg/kk
minimize           1.0e-12 1.0e-12 100000000 100000000
write_data         min.data

#Thermo------------------------------------------------------------------------
thermo             ${thermostep}
thermo_style       custom step temp etotal pe ke vol press density

#Species related---------------------------------------------------------------
fix			   	   rbo all reaxff/bonds/kk ${dumpstep} bonds.reaxc
fix   		   	   spe all reaxff/species/kk 1 1 ${dumpstep} species.out element H C O

#Velocity Command -------------------------------------------------------------
velocity		   all create ${T0} ${seed} dist gaussian

#NVT -------for energy equilibration-------------------------------------------
fix		   		   1 all temp/berendsen ${T0} ${T} $(100.0*dt)
fix 		   	   2 all nve/kk
dump               3 all atom ${dumpstep} ${sim_name}.nvt.lammpstrj
run                ${runstep}
undump             3
unfix              2
unfix	           1
write_restart      ${sim_name}.nvt.restart
write_data         ${sim_name}.nvt.data

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
'''

reinitialize_variables='''
#delete and re-initialize variables-------------------------------------------
variable			runtime delete
variable			runstep delete
variable			resstep	delete
variable			runtime equal {NPT_runtime}
'''.format(**user_inputs)

body_npt = '''
#Re-define variables---------------------------------------------------
variable			runstep equal 1000*${runtime}/${tstep}
variable			resstep equal ${runstep}/${res_freq}

#NPT for density equilibration------------------------------------------------
fix					1 all npt/kk temp ${T0} ${T0} $(100.0*dt) iso ${P0} ${P0} $(1000.0*dt)
dump               	2 all atom ${dumpstep} ${sim_name}.npt.lammpstrj
run                	${runstep}
undump            	2
unfix             	1
write_restart     	${sim_name}.npt.restart
write_data        	${sim_name}.npt.data
'''

submit = '''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ntask}
#SBATCH --time={time}
#SBATCH --mail-user=sahmed73@ucmerced.edu 
#SBATCH --mail-type=ALL
#SBATCH --export=ALL


#-------------------General Preparation-------------------
echo "Nodelist: $SLURM_JOB_NODELIST"
echo "Number of Nodes: $SLURM_NNODES"
echo "Number of Cores per Node: $SLURM_NTASKS_PER_NODE"
module purge
LMP=/data/sahmed73/LAMMPS/LATEST_KOKKOS_PYTHON/lammps-29Sep2021/build/lmp_kk
DIR=input.in
module load openmpi/4.0.6-intel-2021.4.0
#module load openmpi/4.1.1-gcc-8.4.1
#module load mvapich2/2.3.6-intel-2021.4.0
#module load mpich/3.4.2-gcc-8.4.1  
#module load openmpi/4.1.1-intel-2021.4.0
#module load openmpi/4.0.6-intel-2021.4.0

#export OMP_NUM_THREADS=1
NSLOTS=$(($SLURM_NNODES*$SLURM_NTASKS_PER_NODE))
pwd_dir=$(pwd)
OMP_NUM_THREADS=1


#-------------------Run Simulations-------------------

# Run simulation
# mpirun -np $NSLOTS $LMP -log lmp.log -in $DIR
mpirun -np $NSLOTS $LMP -k on -sf kk -in $DIR > output.out
'''.format(**user_inputs)

submit_merced = '''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition_merced}
#SBATCH --nodes={nodes_merced}
#SBATCH --ntasks-per-node={ntask_merced}
#SBATCH --time={time_merced}
#SBATCH --mail-user=sahmed73@ucmerced.edu 
#SBATCH --mail-type=ALL
#SBATCH --export=ALL

#-------------------General Preparation-------------------
echo "Nodelist: $SLURM_JOB_NODELIST"
echo "Number of Nodes: $SLURM_NNODES"
echo "Number of Cores per Node: $SLURM_NTASKS_PER_NODE"

module purge 
module load openmpi/4.0.6-intel-2021.4.0

LMP=/data/sahmed73/LAMMPS/LATEST_KOKKOS_PYTHON/lammps-29Sep2021/build/lmp_kk
DIR=input.in

#module load openmpi/4.1.1-gcc-8.4.1
#module load mvapich2/2.3.6-intel-2021.4.0
#module load mpich/3.4.2-gcc-8.4.1  
#module load openmpi/4.1.1-intel-2021.4.0
#module load openmpi/4.0.6-intel-2021.4.0

#export OMP_NUM_THREADS=1
NSLOTS=$(($SLURM_NNODES*$SLURM_NTASKS_PER_NODE))
pwd_dir=$(pwd)
OMP_NUM_THREADS=1


#-------------------Run Simulations-------------------

# Run simulation
mpirun -np $NSLOTS $LMP -k on -sf kk -in $DIR > output.out
'''.format(**user_inputs)



dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\OPLSAA\PAO_Antioxidant\PAO+A\10PAO+3A\Equilibrate'

sub_dirr = "{NVT_runtime}NVT+{NPT_runtime}NPT".format(**user_inputs)

full_path = os.path.join(dirr, sub_dirr)


if not os.path.exists(full_path):
    os.makedirs(full_path)

# input.in
file_path = os.path.join(full_path, 'input.in')
if os.path.exists(file_path):
    print(f"The file '{file_path}' already exists. Not overwritten.")
else:
    with open(file_path, 'w', newline='\n') as file:
        file.write(initialize_variables)
        file.write(body_nvt)
        file.write(reinitialize_variables)
        file.write(body_npt)

# submit.sh        
file_path = os.path.join(full_path, 'submit.sh')
if os.path.exists(file_path):
    print(f"The file '{file_path}' already exists. Not overwritten.")
else:
    with open(file_path, 'w', newline='\n') as file:
        file.write(submit)

# submit.merced.sh        
file_path = os.path.join(full_path, 'submit.merced.sh')
if os.path.exists(file_path):
    print(f"The file '{file_path}' already exists. Not overwritten.")
else:
    with open(file_path, 'w', newline='\n') as file:
        file.write(submit_merced)
        
print('Successful')