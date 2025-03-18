'''
https://pubs.acs.org/doi/10.1021/jo034808o
'''


from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def calculate_vdw_volume(smiles):
    """
    Calculate the van der Waals volume of a molecule using Bondi's atomic contributions and the Zhao et al. method.
    """
    # Bondi van der Waals volumes (Å³)
    bondi_vdw_volumes = {
        'H': 7.24, 'C': 20.58, 'N': 15.60, 'O': 14.71,
        'F': 13.31, 'Cl': 22.45, 'Br': 26.52, 'I': 32.52,
        'P': 24.43, 'S': 24.43, 'As': 26.52, 'B': 40.48,
        'Si': 38.79, 'Se': 28.73, 'Te': 36.62
    }
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    mol = Chem.AddHs(mol)
    
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    print(atom_counts)
    total_atomic_vdw = sum(bondi_vdw_volumes[atom] * count for atom, count in atom_counts.items() if atom in bondi_vdw_volumes)
    
    num_bonds = mol.GetNumAtoms()-1
    print(num_bonds)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_non_aromatic_rings = num_rings - num_aromatic_rings
    print('# Aromatic ring:',num_aromatic_rings)
    print('# Non-aromatic ring:',num_non_aromatic_rings)
    
    vdw_volume = total_atomic_vdw - (5.92 * num_bonds) - (14.7 * num_aromatic_rings) - (3.8 * num_non_aromatic_rings)
    
    return vdw_volume

# Example Usage
smiles_example = "NC(=O)NN=Cc1ccc(o1)N(=O)=O"
vdw_volume = calculate_vdw_volume(smiles_example)
print(f"Van der Waals Volume: {vdw_volume:.2f} Å³")