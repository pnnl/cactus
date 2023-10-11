from langchain.tools import BaseTool
import pubchempy as pcp


class name_to_SMILES(BaseTool):
    name = "name_to_SMILES"
    description = "Convert the input name into its corresponding SMILES notation"

    def _run(self, input_name: str) -> str:
        """
        Convert chemical name into SMILES notation.

        Parameters:
        input_name (str): The name of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """
        # checking if the input is a string variable or not
        if isinstance(input_name, str):
            # checking for the first compund with the given input name
            result = pcp.get_compounds(str(input_name), "name")[0]
            # print (result)
            return result.isomeric_smiles
        else:
            # if the input is not a string error will be raised
            raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
