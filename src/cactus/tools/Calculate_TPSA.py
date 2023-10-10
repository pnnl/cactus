from rdkit.Chem import MolFromSmiles, Descriptors
from langchain.tools import BaseTool


class calculate_TPSA(BaseTool):
    name = "calculate_TPSA"
    description = "Compute the topological surface area (TPSA) of the given molecule."

    def _run(self, compound: str) -> float:
        """
        Compute the topological surface area (TPSA) of the given molecule.

        Parameters:
        compound: Compound in SMILES format

        Returns:
        float: The TPSA in angstroms^2
        """
        return Descriptors.TPSA(MolFromSmiles(compound))

    async def _arun(self, compound: str) -> float:
        """Use the calculate_TPSA tool asynchronously."""
        raise NotImplementedError()
