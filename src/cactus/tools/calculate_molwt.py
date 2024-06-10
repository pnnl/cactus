"""Tool to calculate the MolWt of a compound."""

from langchain.tools import BaseTool
from rdkit.Chem import Descriptors, MolFromSmiles

DESC = """
Use this tool when you need to calculate the molecular weight of a SMILES string. Units in Dalton.
"""


class CalculateMolWt(BaseTool):
    """Calculate the MolWt of the compound."""

    name = "CalculateMolecularWeight"
    description = DESC

    def _run(self, compound: str) -> float:
        """Compute the exact molecular weight of the given molecule.

        Parameters
        ----------
        compound: Compound in SMILES format

        Returns
        -------
        float: The exact molecular weight in daltons
        """
        mol = MolFromSmiles(compound)
        return Descriptors.ExactMolWt(mol)

    async def _arun(self, compound: str) -> float:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
