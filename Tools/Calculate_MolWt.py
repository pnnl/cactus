#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:10:38 2023

@author: sank064
"""
from rdkit.Chem import MolFromSmiles, Descriptors
from langchain.tools import BaseTool

class calculate_MolWt(BaseTool):
 
    def _run(compound: str) -> float:
        """
        Compute the exact molecular weight of the given molecule.

        Parameters:
        compound: Compound in SMILES format

        Returns:
        float: The exact molecular weight in daltons
        """
        return Descriptors.ExactMolWt(MolFromSmiles(compound))
    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()   
