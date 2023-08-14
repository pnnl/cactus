#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:21:35 2023

@author: sank064
"""
from langchain.tools import BaseTool
from rdkit import Chem
from adme_pred import ADME

class brenk_filter(BaseTool):
    def _run(compound_smiles: str) -> bool:
        """
        Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master

            From the adme-pred-py documentation:
        Brenk (2008) Lessons Learnt from Assembling Screening Libraries for
            Drug Discovery for Neglected Diseases

            Brenk's Structural Alert filter finds fragments "putatively toxic,
            chemically reactive, metabolically unstable or to bear properties
            responsible for poor pharmacokinetics."

            Returns:
                Boolean of whether the molecule triggers the Brenk filter.
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        return mol.brenk()
    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()

#example use
sample = brenk_filter("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
print(sample)