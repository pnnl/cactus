<<<<<<< HEAD
=======

>>>>>>> 95cc70a (refactoring to src/cactus)
from langchain.tools import BaseTool
import pubchempy as pcp


<<<<<<< HEAD
class molecular_formula_to_SMILES(BaseTool):
    name = "molecular_formula_to_SMILES"
    description = "Convert the input molecular formula into its corresponding SMILES notation"

    def _run(self, input_formula: str) -> str:
=======
class molecular_formula_to_SMILES(BaseTool): 
    name="molecular_formula_to_SMILES"
    description="Convert the input molecular formula into its corresponding SMILES notation"
    
    
    def _run(
        self,input_formula: str)-> str:
>>>>>>> 95cc70a (refactoring to src/cactus)
        """
        Convert molecular formula into SMILES notation.

        Parameters:
        input_formula (str): The molecular formula of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """
<<<<<<< HEAD
        # checking if the input is a string variable or not
        if isinstance(input_formula, str):
            # checking for the first compund with the given molecular formula
            result = pcp.get_compounds(str(input_formula), "formula")[0]
            # print (result)
            return result.isomeric_smiles
        else:
            # if the input is not a string error will be raised
=======
    #checking if the input is a string variable or not
        if(isinstance(input_formula,str)):
    #checking for the first compund with the given molecular formula
            result=pcp.get_compounds(str(input_formula),'formula')[0]
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
