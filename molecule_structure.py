# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Wed Oct 25 15:24:54 2023
"""

from rdkit import Chem
from rdkit.Chem import Draw

input_smiles = 'CC'
mol = Chem.MolFromSmiles(input_smiles)
Draw.MolToImage(mol)