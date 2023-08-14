#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:18:18 2023

@author: sank064
"""

from rdkit.Chem import MolFromSmiles, Descriptors
from langchain.tools import BaseTool

class calculate_QED(BaseTool):
    def _run(compound: str) -> float:
        """
        Compute Quantitative Estimate of Druglikeness (QED) of the given molecule. Bickerton et al 2012.

        Parameters:
        compound: Compound in SMILES format

        Returns:
        float: The QED from 0 (druglike) to 1 (not)
        """
        return Descriptors.qed(MolFromSmiles(compound))
    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()    