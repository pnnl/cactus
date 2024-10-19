"""Tool to calculate the Blood Brain Barrier Permeability of a compound."""

from adme_py import ADME
from langchain.tools import BaseTool
from rdkit import Chem


class CalculateBBBPermeant(BaseTool):
    """Tool to calculate the BBB Permeability."""

    name: str = "CalculateBBBPermeant"
    description: str = "calculates the Blood Brain Barrier Permeability of the compound"

    def _run(self, compound_smiles: str) -> str:
        """Calculate Blood Brain Barrier Permeability.

        Returns
        -------
            Boolean of whether the molecule is blood brain permeant or not.

        Notes
        -----
        The BOILED-Egg model is described in:
        A. Daina and V. Zoete, ChemMedChem 2016, 11, 1117.
        https:/doi.org/10.1002/cmdc.201600182
        """
        summary = ADME(compound_smiles).calculate()
        return summary["pharmacokinetics"]["blood_brain_barrier_permeant"]

    async def _arun(self, compound_smiles: str) -> str:
        """Use the calculate_bbb_permeant tool asynchronously."""
        raise NotImplementedError()
