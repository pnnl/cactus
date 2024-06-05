"""Tool to convert a chemical formula into SMILES."""

import pubchempy as pcp
from langchain.tools import BaseTool

DESC = """
Use this tool when you need to convert a chemical formula into a SMILES string.

A chemical formula is a sequence of atomic symbols and/or numeric subscripts that represents the
composition of a compound.

An example of a chemical formula is C6H12O6 or CH4

Only accept input that appears in this format.
"""


class MolecularFormulaToSMILES(BaseTool):
    """Convert Molecular Formula to SMILES."""

    name = "MolecularFormulaToSMILES"
    description = DESC

    def _run(self, input_formula: str) -> str:
        """Convert molecular formula into SMILES notation.

        Parameters
        ----------
        input_formula (str): The molecular formula of the chemical compound.

        Returns
        -------
        str: The SMILES notation in the output format.
        """
        # checking if the input is a string variable or not
        if isinstance(input_formula, str):
            # checking for the first compund with the given molecular formula
            result = pcp.get_compounds(str(input_formula), "formula")[0]

            return result.canonical_smiles

        # if the input is not a string error will be raised
        error_message = f"Invalid input: {input_formula!r}"
        raise ValueError(error_message)

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
