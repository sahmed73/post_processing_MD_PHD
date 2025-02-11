# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Nov 23 03:48:01 2023
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import time

mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)

mb = Chem.MolToMolBlock(mol)
view = py3Dmol.view(width=400, height=400)
view.addModel(mb, 'mol')
view.setStyle({'stick': {}})
view.zoomTo()

# Display the viewer
view.show()
png = view.png()

import base64
with open("molecule_image.png", "wb") as file:
    file.write(base64.b64decode(png.split(',')[1]))