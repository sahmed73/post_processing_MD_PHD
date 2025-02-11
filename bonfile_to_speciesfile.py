# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Dec  8 11:53:31 2023
"""
import magnolia.MD_Converter as conv
import numpy as np
import pandas as pd
import magnolia.bondfile_parser as bfp
import magnolia.speciesfile_parser as sfp


base_oil = 'PAO4'
cutoff=0.5
sim = 'Sim-1'
for cutoff in np.array(range(35,81,5))/100:
    for sim in ['Sim-1','Sim-1','Sim-3']:
        directory = r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Base_Oil\{}\25_{}_200_O2_Soria\Production\1600\{}".format(base_oil,base_oil,sim)
        
        bondfile = directory+'\\bonds.reaxc'
        atomsymbols = ['H','C','O']   
        
        
        df = bfp.get_species_count(bondfile,atomsymbols,cutoff=cutoff)
        speciesfile = bondfile[:bondfile.rfind('\\')]+f'\\species_{cutoff}.out'
        with open(speciesfile,'w') as sf:
            for header in df.columns:
                column = df.loc[:,header]
                # Count the number of species with non-zero values
                column_as_int = column.round(0).astype(int)
                non_zero_species = column_as_int[column_as_int > 0]
                no_moles = non_zero_species.sum()
                no_species = non_zero_species.count()        
                
                # Writing the header line
                string = "# Timestep\tNo_Moles\tNo_Specs\t"
                line = string+"\t".join(non_zero_species.index) + "\n"
                # print(line,end='')
                sf.write(line)
                # Writing the species counts line
                
                species_counts = "\t".join(map(str, non_zero_species.values))
                line = f"{header}\t{no_moles}\t{no_species}\t{species_counts}\n"
                # print(line,end='')
                sf.write(line)
