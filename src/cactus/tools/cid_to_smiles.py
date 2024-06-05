"""Tool to calculate CID to SMILES."""

import pubchempy as pcp
from langchain.tools import BaseTool


class CidToSMILES(BaseTool):
    """Convert Cid to SMILES."""

    name = "CidToSMILES"
    description = "Convert the input Pubchem id into its corresponding SMILES notation"

    def _run(self, input_id: str) -> str:
        """Convert PubChem id into SMILES notation.

        Parameters
        ----------
        input_id (str): The PubChem id of the chemical compound.

        Returns
        -------
        str: The SMILES notation in the output format.
        """
        # checking if the input is a integer variable or not
        if isinstance(input_id, str):
            result = pcp.Compound.from_cid(int(input_id))
            return result.canonical_smiles

        # if the input is not a string error will be raised
        error_message = f"Invalid input: {input_id!r}"
        raise ValueError(error_message)

    async def _arun(self, input_id: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
