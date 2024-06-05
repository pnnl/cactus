"""Tool for converting common name to smiles."""

import pubchempy as pcp
from langchain.tools import BaseTool

DESC = """
Use this tool when you need to convert a molecule's common name to it's corresponding SMILES string.

A common name could look like this: Acetate, Sulfate, Alkene, Benzene.

Most generic searches will be in this format.
"""


class NameToSMILES(BaseTool):
    """Convert chemical name into SMILES notation."""

    name = "CommonNameToSMILES"
    description = DESC

    def _run(self, input_name: str) -> str:
        """Convert chemical name into SMILES notation.

        Parameters
        ----------
             input_name (str): The name of the chemical compound.

        Returns
        -------
             str: The SMILES notation in the output format.
        """
        # checking if the input is a string variable or not
        if isinstance(input_name, str):
            # checking for the first compund with the given input name
            result = pcp.get_compounds(str(input_name), "name")[0]
            return result.canonical_smiles

        # if the input is not a string error will be raised
        error_message = f"Invalid input: {input_name!r}"
        raise ValueError(error_message)

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
