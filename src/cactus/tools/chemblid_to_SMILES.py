from langchain.tools import BaseTool
from chembl_webresource_client.new_client import new_client


class chemblid_to_SMILES(BaseTool):
    name = "chemblid_to_SMILES"
    description = "Convert the input chemblid into its corresponding SMILES notation"

    def _run(self, input_id: str) -> str:
        """
        Convert ChEMBL id into SMILES notation.

        Parameters:
        input_id (str): The ChEMBL database id of the chemical compound.

        Returns:
        str: The SMILES notation in the output format.
        """
        if isinstance(input_id, str):
            # from chembl_webresource_client.new_client import new_client
            molecule = new_client.molecule
            molecule.set_format("json")
            record_via_client = molecule.get(input_id)
            smiles_from_json = record_via_client["molecule_structures"]["canonical_smiles"]
            return smiles_from_json
        else:
            raise ValueError("Invalid input")

    async def _arun(self, input_name: str) -> str:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
