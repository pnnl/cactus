from langchain.tools import BaseTool
import pubchempy as pcp


class inchikey_to_SMILES(BaseTool):
    name = "convert_to_SMILES"
    description = "Convert the input inchikey into its corresponding SMILES notation"

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
            # print (result)
            return result.isomeric_smiles
        else:
            raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
