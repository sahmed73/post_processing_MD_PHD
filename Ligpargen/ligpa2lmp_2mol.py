import numpy as np
import pandas as pd

# Define helper functions to import data from files
def num_and_types(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0], header=None).to_numpy().flatten()

def import_dims(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1], header=None).to_numpy()

def import_masses(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1], header=None).to_numpy()

def import_pair_coeffs(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1,2], header=None).to_numpy()

def import_bond_coeffs(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1,2], header=None).to_numpy()

def import_angle_coeffs(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1,2], header=None).to_numpy()

def import_dihedral_coeffs(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1,2,3,4], header=None).to_numpy()

def import_improper_coeffs(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[0,1,2,3], header=None).to_numpy()

def import_atoms(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, header=None).to_numpy()

def import_bonds(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, header=None).to_numpy()

def import_angles(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, header=None).to_numpy()

def import_torsional(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, header=None).to_numpy()

def import_packmol(filename, data_lines):
    return pd.read_csv(filename, sep=" ", skiprows=data_lines[0]-1, nrows=data_lines[1]-data_lines[0]+1, usecols=[6,7,8], header=None).to_numpy()

# Assign file paths and number of molecules
filename_a1 = "C:/Users/arup2/OneDrive - University of California Merced/Desktop/LAMMPS/OPLSAA/PAO_Antioxidant/PAO+B/10PAO+3B/DataFile/PAO.lmp"
filename_a2 = "C:/Users/arup2/OneDrive - University of California Merced/Desktop/LAMMPS/OPLSAA/PAO_Antioxidant/PAO+B/10PAO+3B/DataFile/B.lmp"
nMols_a1 = 10
nMols_a2 = 3

sq = 28.000
boxsize = np.array([[0.0, 0.0, 0.0], [sq, sq, sq]])
packed_pdbfile = "C:/Users/arup2/OneDrive - University of California Merced/Desktop/LAMMPS/OPLSAA/PAO_Antioxidant/PAO+B/10PAO+3B/DataFile/10PAO_3B.pdb"

# Description of molecules
nums_a1 = num_and_types(filename_a1, [3, 7])
types_a1 = num_and_types(filename_a1, [9, 13])
dims_a1 = import_dims(filename_a1, [15, 17])

nums_a2 = num_and_types(filename_a2, [3, 7])
types_a2 = num_and_types(filename_a2, [9, 13])
dims_a2 = import_dims(filename_a2, [15, 17])

nsections = 14

# Initialize data lines
data_lines_a1 = np.zeros((nsections, 2), dtype=int)
data_lines_a2 = np.zeros((nsections, 2), dtype=int)

# Populate data lines
for i in range(1, nsections + 1):
    if i == 1:
        data_lines_a1[i-1] = [3, 3 + len(nums_a1) - 1]
        data_lines_a2[i-1] = [3, 3 + len(nums_a2) - 1]
    elif i == 2:
        data_lines_a1[i-1] = [data_lines_a1[i-2, 1] + 2, data_lines_a1[i-2, 1] + 2 + len(nums_a1) - 1]
        data_lines_a2[i-1] = [data_lines_a2[i-2, 1] + 2, data_lines_a2[i-2, 1] + 2 + len(nums_a2) - 1]
    elif i == 3:
        data_lines_a1[i-1] = [data_lines_a1[i-2, 1] + 2, data_lines_a1[i-2, 1] + 2 + len(dims_a1) - 1]
        data_lines_a2[i-1] = [data_lines_a2[i-2, 1] + 2, data_lines_a2[i-2, 1] + 2 + len(dims_a2) - 1]
    else:
        emptyrow = 3 if i > 3 else 1
        if i in [4, 5, 6, 7, 8, 9]:
            num_types = types_a1[i-4] if i < 10 else nums_a1[i-10]
        else:
            num_types = nums_a1[i-10] if i < 15 else nums_a1[i-15]
        
        if np.isnan(num_types):
            print(f"Error: num_types is NaN at section {i}")
            continue
        
        data_lines_a1[i-1] = [data_lines_a1[i-2, 1] + emptyrow + 1, data_lines_a1[i-2, 1] + emptyrow + num_types]
        num_types = types_a2[i-4] if i < 10 else nums_a2[i-10]
        data_lines_a2[i-1] = [data_lines_a2[i-2, 1] + emptyrow + 1, data_lines_a2[i-2, 1] + emptyrow + num_types]

# Check for NaNs in data lines
if np.any(np.isnan(data_lines_a1)) or np.any(np.isnan(data_lines_a2)):
    print("NaN detected in data lines. Exiting.")
else:
    # Import coefficients
    masses_a1 = import_masses(filename_a1, data_lines_a1[3])
    pair_coeffs_a1 = import_pair_coeffs(filename_a1, data_lines_a1[4])
    bond_coeffs_a1 = import_bond_coeffs(filename_a1, data_lines_a1[5])
    angle_coeffs_a1 = import_angle_coeffs(filename_a1, data_lines_a1[6])
    dihedral_coeffs_a1 = import_dihedral_coeffs(filename_a1, data_lines_a1[7])
    improper_coeffs_a1 = import_improper_coeffs(filename_a1, data_lines_a1[8])

    masses_a2 = import_masses(filename_a2, data_lines_a2[3])
    pair_coeffs_a2 = import_pair_coeffs(filename_a2, data_lines_a2[4])
    bond_coeffs_a2 = import_bond_coeffs(filename_a2, data_lines_a2[5])
    angle_coeffs_a2 = import_angle_coeffs(filename_a2, data_lines_a2[6])
    dihedral_coeffs_a2 = import_dihedral_coeffs(filename_a2, data_lines_a2[7])
    improper_coeffs_a2 = import_improper_coeffs(filename_a2, data_lines_a2[8])

    # Import atoms, bonds, angles, dihedrals, and impropers sections
    atoms_a1 = import_atoms(filename_a1, data_lines_a1[9])
    bonds_a1 = import_bonds(filename_a1, data_lines_a1[10])
    angles_a1 = import_angles(filename_a1, data_lines_a1[11])
    dihedrals_a1 = import_torsional(filename_a1, data_lines_a1[12])
    impropers_a1 = import_torsional(filename_a1, data_lines_a1[13])

    atoms_a2 = import_atoms(filename_a2, data_lines_a2[9])
    bonds_a2 = import_bonds(filename_a2, data_lines_a2[10])
    angles_a2 = import_angles(filename_a2, data_lines_a2[11])
    dihedrals_a2 = import_torsional(filename_a2, data_lines_a2[12])
    impropers_a2 = import_torsional(filename_a2, data_lines_a2[13])

    # Import coordinates from packmol file
    atomsXYZ = import_packmol(packed_pdbfile, [6, 6 + nums_a1[0]*nMols_a1 + nums_a2[0]*nMols_a2 - 1])

    # Identify the Identical coefficients and eliminate them
    def unique_coeffs(coeffs):
        unique_rows, indices, inverse_indices = np.unique(coeffs[:, 1:], axis=0, return_index=True, return_inverse=True)
        new_coeffs = np.hstack((np.arange(1, len(unique_rows) + 1)[:, np.newaxis], unique_rows))
        return new_coeffs, indices, inverse_indices + 1

    pair_coeffs2_a1, ia_a1, ic_a1 = unique_coeffs(pair_coeffs_a1)
    atoms_a1[:, 2] = ic_a1
    masses2_a1 = np.hstack((np.arange(1, len(pair_coeffs2_a1) + 1)[:, np.newaxis], masses_a1[ia_a1, 1:]))

    bond_coeffs2_a1, ia_a1, ic_a1 = unique_coeffs(bond_coeffs_a1)
    bonds_a1[:, 1] = ic_a1

    angle_coeffs2_a1, ia_a1, ic_a1 = unique_coeffs(angle_coeffs_a1)
    angles_a1[:, 1] = ic_a1

    dihedral_coeffs2_a1, ia_a1, ic_a1 = unique_coeffs(dihedral_coeffs_a1)
    dihedrals_a1[:, 1] = ic_a1

    improper_coeffs2_a1, ia_a1, ic_a1 = unique_coeffs(improper_coeffs_a1)
    impropers_a1[:, 1] = ic_a1

    pair_coeffs2_a2, ia_a2, ic_a2 = unique_coeffs(pair_coeffs_a2)
    atoms_a2[:, 0] += nMols_a1 * max(atoms_a1[:, 0])
    atoms_a2[:, 1] += nMols_a1
    atoms_a2[:, 2] = ic_a2 + max(atoms_a1[:, 2])
    pair_coeffs2_a2[:, 0] += max(pair_coeffs2_a1[:, 0])
    masses2_a2 = np.hstack((np.arange(1, len(pair_coeffs2_a2) + 1)[:, np.newaxis] + max(masses2_a1[:, 0]), masses_a2[ia_a2, 1:]))

    bond_coeffs2_a2, ia_a2, ic_a2 = unique_coeffs(bond_coeffs_a2)
    bonds_a2[:, 0] += nMols_a1 * max(bonds_a1[:, 0])
    bonds_a2[:, 1] = ic_a2 + max(bonds_a1[:, 1])
    bonds_a2[:, 2:] += nMols_a1 * max(atoms_a1[:, 0])
    bond_coeffs2_a2[:, 0] += max(bond_coeffs2_a1[:, 0])

    angle_coeffs2_a2, ia_a2, ic_a2 = unique_coeffs(angle_coeffs_a2)
    angles_a2[:, 0] += nMols_a1 * max(angles_a1[:, 0])
    angles_a2[:, 1] = ic_a2 + max(angles_a1[:, 1])
    angles_a2[:, 2:] += nMols_a1 * max(atoms_a1[:, 0])
    angle_coeffs2_a2[:, 0] += max(angle_coeffs2_a1[:, 0])

    dihedral_coeffs2_a2, ia_a2, ic_a2 = unique_coeffs(dihedral_coeffs_a2)
    dihedrals_a2[:, 0] += nMols_a1 * max(dihedrals_a1[:, 0])
    dihedrals_a2[:, 1] = ic_a2 + max(dihedrals_a1[:, 1])
    dihedrals_a2[:, 2:] += nMols_a1 * max(atoms_a1[:, 0])
    dihedral_coeffs2_a2[:, 0] += max(dihedral_coeffs2_a1[:, 0])

    improper_coeffs2_a2, ia_a2, ic_a2 = unique_coeffs(improper_coeffs_a2)
    impropers_a2[:, 0] += nMols_a1 * max(impropers_a1[:, 0])
    impropers_a2[:, 1] = ic_a2 + max(impropers_a1[:, 1])
    impropers_a2[:, 2:] += nMols_a1 * max(atoms_a1[:, 0])
    improper_coeffs2_a2[:, 0] += max(improper_coeffs2_a1[:, 0])

    # Develop atoms, bonds, angles, dihedrals, and impropers for nMols of molecules
    total_atoms = nums_a1[0]*nMols_a1 + nums_a2[0]*nMols_a2
    total_bonds = nums_a1[1]*nMols_a1 + nums_a2[1]*nMols_a2
    total_angles = nums_a1[2]*nMols_a1 + nums_a2[2]*nMols_a2
    total_dihedrals = nums_a1[3]*nMols_a1 + nums_a2[3]*nMols_a2
    total_impropers = nums_a1[4]*nMols_a1 + nums_a2[4]*nMols_a2

    atoms2_all = np.zeros((total_atoms, atoms_a2.shape[1]))
    bonds2_all = np.zeros((total_bonds, bonds_a2.shape[1]))
    angles2_all = np.zeros((total_angles, angles_a2.shape[1]))
    dihedrals2_all = np.zeros((total_dihedrals, dihedrals_a2.shape[1]))
    impropers2_all = np.zeros((total_impropers, impropers_a2.shape[1]))

    atoms2_all[:, 0] = np.arange(1, total_atoms + 1)
    bonds2_all[:, 0] = np.arange(1, total_bonds + 1)
    angles2_all[:, 0] = np.arange(1, total_angles + 1)
    dihedrals2_all[:, 0] = np.arange(1, total_dihedrals + 1)
    impropers2_all[:, 0] = np.arange(1, total_impropers + 1)

    for i in range(nMols_a1):
        start_row_a1 = i * nums_a1
        end_row_a1 = (i + 1) * nums_a1
        if i == 0:
            atoms2_all[start_row_a1[0]:end_row_a1[0], :] = atoms_a1
            bonds2_all[start_row_a1[1]:end_row_a1[1], :] = bonds_a1
            angles2_all[start_row_a1[2]:end_row_a1[2], :] = angles_a1
            dihedrals2_all[start_row_a1[3]:end_row_a1[3], :] = dihedrals_a1
            impropers2_all[start_row_a1[4]:end_row_a1[4], :] = impropers_a1
        else:
            atoms2_all[start_row_a1[0]:end_row_a1[0], 1:] = [atoms2_all[start_row_a1[0]-nums_a1[0]:end_row_a1[0]-nums_a1[0], 1]+1, atoms2_all[start_row_a1[0]-nums_a1[0]:end_row_a1[0]-nums_a1[0], 2:]]
            bonds2_all[start_row_a1[1]:end_row_a1[1], 1:] = [bonds2_all[start_row_a1[1]-nums_a1[1]:end_row_a1[1]-nums_a1[1], 1], bonds2_all[start_row_a1[1]-nums_a1[1]:end_row_a1[1]-nums_a1[1], 2:]+nums_a1[0]]
            angles2_all[start_row_a1[2]:end_row_a1[2], 1:] = [angles2_all[start_row_a1[2]-nums_a1[2]:end_row_a1[2]-nums_a1[2], 1], angles2_all[start_row_a1[2]-nums_a1[2]:end_row_a1[2]-nums_a1[2], 2:]+nums_a1[0]]
            dihedrals2_all[start_row_a1[3]:end_row_a1[3], 1:] = [dihedrals2_all[start_row_a1[3]-nums_a1[3]:end_row_a1[3]-nums_a1[3], 1], dihedrals2_all[start_row_a1[3]-nums_a1[3]:end_row_a1[3]-nums_a1[3], 2:]+nums_a1[0]]
            impropers2_all[start_row_a1[4]:end_row_a1[4], 1:] = [impropers2_all[start_row_a1[4]-nums_a1[4]:end_row_a1[4]-nums_a1[4], 1], impropers2_all[start_row_a1[4]-nums_a1[4]:end_row_a1[4]-nums_a1[4], 2:]+nums_a1[0]]

    for i in range(nMols_a2):
        start_row_a2 = nums_a1[0]*nMols_a1 + i * nums_a2
        end_row_a2 = nums_a1[0]*nMols_a1 + (i + 1) * nums_a2
        if i == 0:
            atoms2_all[start_row_a2[0]:end_row_a2[0], :] = atoms_a2
            bonds2_all[start_row_a2[1]:end_row_a2[1], :] = bonds_a2
            angles2_all[start_row_a2[2]:end_row_a2[2], :] = angles_a2
            dihedrals2_all[start_row_a2[3]:end_row_a2[3], :] = dihedrals_a2
            impropers2_all[start_row_a2[4]:end_row_a2[4], :] = impropers_a2
        else:
            atoms2_all[start_row_a2[0]:end_row_a2[0], 1:] = [atoms2_all[start_row_a2[0]-nums_a2[0]:end_row_a2[0]-nums_a2[0], 1]+1, atoms2_all[start_row_a2[0]-nums_a2[0]:end_row_a2[0]-nums_a2[0], 2:]]
            bonds2_all[start_row_a2[1]:end_row_a2[1], 1:] = [bonds2_all[start_row_a2[1]-nums_a2[1]:end_row_a2[1]-nums_a2[1], 1], bonds2_all[start_row_a2[1]-nums_a2[1]:end_row_a2[1]-nums_a2[1], 2:]+nums_a2[0]]
            angles2_all[start_row_a2[2]:end_row_a2[2], 1:] = [angles2_all[start_row_a2[2]-nums_a2[2]:end_row_a2[2]-nums_a2[2], 1], angles2_all[start_row_a2[2]-nums_a2[2]:end_row_a2[2]-nums_a2[2], 2:]+nums_a2[0]]
            dihedrals2_all[start_row_a2[3]:end_row_a2[3], 1:] = [dihedrals2_all[start_row_a2[3]-nums_a2[3]:end_row_a2[3]-nums_a2[3], 1], dihedrals2_all[start_row_a2[3]-nums_a2[3]:end_row_a2[3]-nums_a2[3], 2:]+nums_a2[0]]
            impropers2_all[start_row_a2[4]:end_row_a2[4], 1:] = [impropers2_all[start_row_a2[4]-nums_a2[4]:end_row_a2[4]-nums_a2[4], 1], impropers2_all[start_row_a2[4]-nums_a2[4]:end_row_a2[4]-nums_a2[4], 2:]+nums_a2[0]]

    atoms2_all[:, 4:] = atomsXYZ

    # Write LAMMPS data file
    output_file = "C:/path/to/output_file.txt"
    with open(output_file, 'w') as file:
        file.write("LAMMPS data file of the system for full atom style\n")
        file.write(f"\n{total_atoms} atoms\n")
        file.write(f"{total_bonds} bonds\n")
        file.write(f"{total_angles} angles\n")
        file.write(f"{total_dihedrals} dihedrals\n")
        file.write(f"{total_impropers} impropers\n")
        file.write(f"\n{len(pair_coeffs2_a1) + len(pair_coeffs2_a2)} atom types\n")
        file.write(f"{len(bond_coeffs2_a1) + len(bond_coeffs2_a2)} bond types\n")
        file.write(f"{len(angle_coeffs2_a1) + len(angle_coeffs2_a2)} angle types\n")
        file.write(f"{len(dihedral_coeffs2_a1) + len(dihedral_coeffs2_a2)} dihedral types\n")
        file.write(f"{len(improper_coeffs2_a1) + len(improper_coeffs2_a2)} improper types\n")
        file.write(f"\n{boxsize[0,0]:.6f} {boxsize[0,1]:.6f} xlo xhi\n")
        file.write(f"{boxsize[1,0]:.6f} {boxsize[1,1]:.6f} ylo yhi\n")
        file.write(f"{boxsize[2,0]:.6f} {boxsize[2,1]:.6f} zlo zhi\n")
        file.write("\n Masses \n\n")
        for mass in np.vstack((masses2_a1, masses2_a2)):
            file.write(f"{mass[0]} {mass[1]:.4f}\n")
        file.write("\n Pair Coeffs \n\n")
        for coeff in np.vstack((pair_coeffs2_a1, pair_coeffs2_a2)):
            file.write(f"{coeff[0]} {coeff[1]:.4f} {coeff[2]:.4f}\n")
        file.write("\n Bond Coeffs \n\n")
        for coeff in np.vstack((bond_coeffs2_a1, bond_coeffs2_a2)):
            file.write(f"{coeff[0]} {coeff[1]:.4f} {coeff[2]:.4f}\n")
        file.write("\n Angle Coeffs \n\n")
        for coeff in np.vstack((angle_coeffs2_a1, angle_coeffs2_a2)):
            file.write(f"{coeff[0]} {coeff[1]:.4f} {coeff[2]:.4f}\n")
        file.write("\n Dihedral Coeffs \n\n")
        for coeff in np.vstack((dihedral_coeffs2_a1, dihedral_coeffs2_a2)):
            file.write(f"{coeff[0]} {coeff[1]:.4f} {coeff[2]:.4f} {coeff[3]:.4f} {coeff[4]:.4f}\n")
        file.write("\n Improper Coeffs \n\n")
        for coeff in np.vstack((improper_coeffs2_a1, improper_coeffs2_a2)):
            file.write(f"{coeff[0]} {coeff[1]:.4f} {coeff[2]:.4f} {coeff[3]:.4f}\n")
        file.write("\n Atoms \n\n")
        for atom in atoms2_all:
            file.write(f"{int(atom[0])} {int(atom[1])} {int(atom[2])} {atom[3]:.5f} {atom[4]:.5f} {atom[5]:.5f} {atom[6]:.5f}\n")
        file.write("\n Bonds \n\n")
        for bond in bonds2_all:
            file.write(f"{int(bond[0])} {int(bond[1])} {int(bond[2])} {int(bond[3])}\n")
        file.write("\n Angles \n\n")
        for angle in angles2_all:
            file.write(f"{int(angle[0])} {int(angle[1])} {int(angle[2])} {int(angle[3])} {int(angle[4])}\n")
        file.write("\n Dihedrals \n\n")
        for dihedral in dihedrals2_all:
            file.write(f"{int(dihedral[0])} {int(dihedral[1])} {int(dihedral[2])} {int(dihedral[3])} {int(dihedral[4])} {int(dihedral[5])}\n")
        file.write("\n Impropers \n\n")
        for improper in impropers2_all:
            file.write(f"{int(improper[0])} {int(improper[1])} {int(improper[2])} {int(improper[3])} {int(improper[4])} {int(improper[5])}\n")

    print('Success!!!')
