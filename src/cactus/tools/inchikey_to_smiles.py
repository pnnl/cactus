"""Tool to convert InchiKey to SMILES."""

import pubchempy as pcp
from langchain.tools import BaseTool

DESC = """
Use this tool when you need to convert a molecule's inchikey
to it's corresponding SMILES string.

Examples of InchiKeys:
- HUMNYLRZRPPJDN-UHFFFAOYSA-N
- BJKAKPMCBYMRJI-UHFFFAOYSA-N
- ZRSNZINYAWTAHE-UHFFFAOYSA-N
"""


class InchikeyToSMILES(BaseTool):
    """Convert Inchikey to SMILES."""

    name = "InchikeyToSMILES"
    description = DESC

    def _run(self, input_inchikey: str) -> str:
        """Convert inchikey into SMILES notation.

        Parameters
        ----------
        input_name (str): The InchIkey of the chemical compound.

        Returns
        -------
        str: The SMILES notation in the output format.
        """
        if isinstance(input_inchikey, str):
            result = pcp.get_compounds(str(input_inchikey), "inchikey")[0]
            return result.canonical_smiles

        error_message = f"Invalid input: {input_inchikey!r}"
        raise ValueError(error_message)

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
