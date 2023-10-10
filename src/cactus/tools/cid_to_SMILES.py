<<<<<<< HEAD
=======
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:23:13 2023

@author: sank064
"""
>>>>>>> 95cc70a (refactoring to src/cactus)
from langchain.tools import BaseTool
import pubchempy as pcp


<<<<<<< HEAD
class cid_to_SMILES(BaseTool):
    name = "cid_to_SMILES"
    description = "Convert the input Pubchem id into its corresponding SMILES notation"

    def _run(self, input_id: int) -> str:
=======
class cid_to_SMILES(BaseTool):     
    name="cid_to_SMILES"
    description="Convert the input Pubchem id into its corresponding SMILES notation"
    
    def _run(
            self,input_id: int)-> str:
>>>>>>> 95cc70a (refactoring to src/cactus)
        """
        Convert PubChem id into SMILES notation.

        Parameters:
        input_id (str): The PubChem id of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """
<<<<<<< HEAD
        # checking if the input is a integer variable or not
        if isinstance(input_id, int):
            result = pcp.Compound.from_cid(int(input_id))
            # print (result)
            return result.isomeric_smiles
        else:
            # if the input is not a string error will be raised
=======
    #checking if the input is a integer variable or not
        if(isinstance(input_id,int)):
            result=pcp.Compound.from_cid(int(input_id))
            #print (result)
            return result.isomeric_smiles
        else:
    #if the input is not a string error will be raised
>>>>>>> 95cc70a (refactoring to src/cactus)
            raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
<<<<<<< HEAD
        raise NotImplementedError()
=======
        raise NotImplementedError()
>>>>>>> 95cc70a (refactoring to src/cactus)
