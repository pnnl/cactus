from langchain.tools import BaseTool
import pubchempy as pcp


class cid_to_SMILES(BaseTool):
    name = "cid_to_SMILES"
    description = "Convert the input Pubchem id into its corresponding SMILES notation"

    def _run(self, input_id: int) -> str:
        """
        Convert PubChem id into SMILES notation.

        Parameters:
        input_id (str): The PubChem id of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """

        # checking if the input is a integer variable or not
        if isinstance(input_id, int):
            result = pcp.Compound.from_cid(int(input_id))
            return result.canonical_smiles

        # if the input is not a string error will be raised
        raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
