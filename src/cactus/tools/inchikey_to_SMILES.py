from langchain.tools import BaseTool
import pubchempy as pcp

DESC = """
Use this tool when you need to convert a molecule's inchikey
to it's corresponding SMILES string.

An InChIKey currently consists of three parts separated by hyphens
, of 14, 10 and one character(s), respectively,
like XXXXXXXXXXXXXX-YYYYYYYYFV-P

Only accept input that appears in this format.
"""


class inchikey_to_SMILES(BaseTool):
    name = "inchikey_to_SMILES"
    description = DESC

    def _run(self, input_inchikey: str) -> str:
        """
        Convert inchikey into SMILES notation.

        Parameters:
        input_name (str): The InchIkey of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """
        if isinstance(input_inchikey, str):
            result = pcp.get_compounds(str(input_inchikey), "inchikey")[0]
            return result.canonical_smiles

        raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
