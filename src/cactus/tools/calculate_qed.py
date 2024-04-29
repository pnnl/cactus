from langchain.tools import BaseTool
from rdkit.Chem import Descriptors, MolFromSmiles


class CalculateQED(BaseTool):
    name = "CalculateQED"
    description = "Compute Quantitative Estimate of Druglikeness (QED) of the given molecule"

    def _run(self, compound: str) -> float:
        """
        Compute Quantitative Estimate of Druglikeness (QED) of the given molecule.
        Bickerton et al 2012.

        Parameters:
        compound: Compound in SMILES format

        Returns:
        float: The QED from 0 (druglike) to 1 (not)
        """
        return Descriptors.qed(MolFromSmiles(compound))

    async def _arun(self, compound: str) -> float:
        """Use the calculate_QED tool asynchronously."""
        raise NotImplementedError()
