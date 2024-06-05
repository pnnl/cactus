"""Tool to calculate the TPSA of a compound."""

from langchain.tools import BaseTool
from rdkit.Chem import MolFromSmiles, rdMolDescriptors


class CalculateTPSA(BaseTool):
    """Calculate the TPSA of the compound."""

    name = "CalculateTPSA"
    description = "Compute the Topological polar surface area (TPSA) of the given molecule."

    def _run(self, compound: str) -> float:
        """Compute the topological surface area (TPSA) of the given molecule.

        Parameters
        ----------
        compound: Compound in SMILES format

        Returns
        -------
        float: The TPSA in angstroms^2
        """
        return rdMolDescriptors.CalcTPSA(MolFromSmiles(compound))

    async def _arun(self, compound: str) -> float:
        """Use the calculate_TPSA tool asynchronously."""
        raise NotImplementedError()
