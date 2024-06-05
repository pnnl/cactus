"""Tool to convert ZincID to SMILES."""

import smilite
from langchain.tools import BaseTool


class ZincIDToSMILES(BaseTool):
    """Convert ZincID to SMILES."""

    name = "zinc_id_to_SMILES"
    description = "Convert the input zinc id into its corresponding SMILES notation"

    def _run(self, input_id: str) -> str:
        """Convert zinc id into SMILES notation.

        Parameters
        ----------
        input_id (str): The ZINC15 id of the chemical compound.

        Returns
        -------
        str: The SMILES notation in the output format.
        """
        # checking if the input is a string variable or not

        if isinstance(input_id, str):
            # checking for the first compund with the given molecular formula
            smile_str = smilite.get_zinc_smile(input_id, backend="zinc15")

            return smile_str

        # if the input is not a string error will be raised
        error_message = f"Invalid input: {input_id!r}"
        raise ValueError(error_message)

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
