
from rdkit.Chem import MolFromSmiles, Descriptors
from langchain.tools import BaseTool

class calculate_MolWt(BaseTool):
    name="calculate_MolWt"
    description="Compute the exact molecular weight of the given molecule." 

    def _run(compound: str) -> float:
        """
        Compute the exact molecular weight of the given molecule.

        Parameters:
        compound: Compound in SMILES format

        Returns:
        float: The exact molecular weight in daltons
        """
        return Descriptors.ExactMolWt(MolFromSmiles(compound))
    
    async def _arun(self, compound: str) -> float:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()   