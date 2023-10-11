#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:23:11 2023

@author: sank064
"""
from langchain.tools import BaseTool
import pubchempy as pcp


class cas_to_SMILES(BaseTool):
    name = "cas_to_SMILES"
    description = "Convert the input cas into its corresponding SMILES notation"

    def _run(self, input_query: str) -> str:
        """
        Convert cas number into SMILES notation.

        Parameters:
        input_query (str): The cas number of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """

        if isinstance(input_query, str):
            result = pcp.get_compounds(str(input_query), "name")[0]
            # print (result)
            return result
        #         get_cid("input_query", from = "xref/RN")
        #         return pc_sect(702,"canonical smiles")
        else:
            raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
