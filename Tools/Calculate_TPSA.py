#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:20:09 2023

@author: sank064
"""
from rdkit.Chem import MolFromSmiles, Descriptors
from langchain.tools import BaseTool

class calculate_TPSA(BaseTool):
    def _run(compound: str) -> float:
        """
        Compute the topological surface area (TPSA) of the given molecule.

        Parameters:
        compound: Compound in SMILES format

        Returns:
        float: The TPSA in angstroms^2
        """
        return Descriptors.TPSA(MolFromSmiles(compound))
    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()