"""Tool to convert CAS to SMILES."""

import pubchempy as pcp
from langchain.tools import BaseTool


class CasToSMILES(BaseTool):
    """Convert Cas to SMILES."""

    name = "CasToSMILES"
    description = "Convert the input cas into its corresponding SMILES notation"

    def _run(self, input_query: str) -> str:
        """Convert cas number into SMILES notation.

        Parameters
        ----------
        input_query (str): The cas number of the chemical compound.

        Returns
        -------
        str: The SMILES notation in the output format.
        """
        if isinstance(input_query, str):
            result = pcp.get_compounds(str(input_query), "name")[0]
            # print (result)
            return result.canonical_smiles

        error_message = f"Invalid input: {input_query!r}"
        raise ValueError(error_message)

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
