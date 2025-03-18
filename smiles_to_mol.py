import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

'''
argv[0] = script name (ignored)
argv[1] = csv_file (input CSV with name and smiles columns)
argv[2] = output_dir (directory where .mol files will be saved)
'''

csv_file = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\Solubility_Automation\test_DELETE\smiles.csv" 
output_dir = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\borgstore\Solubility\Solubility_Automation\test_DELETE"

# Check if correct number of arguments are provided
if len(sys.argv) > 1:
    if len(sys.argv) < 3:
        print("Usage: python script.py <csv_file> <output_dir>")
        
    # Get input arguments
    csv_file = sys.argv[1] 
    output_dir = sys.argv[2]
    
    sys.exit(1)

# if os.path.exists(output_dir):
#     print(f"Error: Output directory '{output_dir}' already exists. Please specify a new directory.")
#     sys.exit(1)
    
# Create output directory if it doesn't exist
# os.makedirs(output_dir)

# Read the CSV file
df = pd.read_csv(csv_file)

# Check if required columns exist
if "name" not in df.columns or "smiles" not in df.columns:
    raise ValueError("CSV must contain 'name' and 'smiles' columns.")

# Loop through each row and generate .mol files
for index, row in df.iterrows():
    name = row["name"]
    smiles = row["smiles"]

    # Convert SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)

    # Skip invalid SMILES
    if mol is None:
        print(f"Warning: Could not parse SMILES for {name}. Skipping...")
        continue

    # Add hydrogens (optional)
    mol = Chem.AddHs(mol)
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    # Define the folder path (parent directory)
    molecule_dir = os.path.join(output_dir, f"{name}", "DataFile")
    os.makedirs(molecule_dir, exist_ok=True)  # Create only directories

    # Construct the .mol file path
    mol_file_path = os.path.join(molecule_dir, f"{name}.mol")

    # Save as MOL file
    Chem.MolToMolFile(mol, mol_file_path)

print("SMILES to .mol conversion completed successfully!")
