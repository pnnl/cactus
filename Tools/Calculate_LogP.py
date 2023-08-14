#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 13:02:13 2023

@author: sank064
"""
from langchain.tools import BaseTool
from rdkit import Chem
from rdkit.Chem import Descriptors

class calculate_logp(BaseTool):
    def _run(compound_smiles: str) -> float:
        """
        RdKit descriptor module to calculate LogP from the atom based Crippen approach
        https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
        https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        return Descriptors.MolLogP(mol)

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError() 
    #example use
sample = calculate_logp("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
print(sample)