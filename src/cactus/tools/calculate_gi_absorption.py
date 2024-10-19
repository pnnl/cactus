"""Tool to calculate the GI Absorption of the compound."""

from adme_py import ADME
from langchain.tools import BaseTool


class CalculateGIAbsorption(BaseTool):
    """Calculate the GI Absorption of a compound."""

    name = "calculate_gi_absorption"
    description = "returns whether the gastrointestinal absorption is high or low"

    def _run(self, compound_smiles: str) -> str:
        """Calculate the gastrointestinal absorption of the molecule.

        Parameters
        ----------
        compound_smiles : str
            The smiles string of the molecule.

        Returns
        -------
            string of whether the gastrointestinal absorption is 'high' or 'low'

        Notes
        -----
        The BOILED-Egg model is described in:
        A. Daina and V. Zoete, ChemMedChem 2016, 11, 1117.
        https:/doi.org/10.1002/cmdc.201600182
        """
        summary = ADME(compound_smiles).calculate()
        return summary["pharmacokinetics"]["gastrointestinal_absorption"]

    async def _arun(self, compound_smiles: str) -> str:
        """Use the calculate_gi_absorption tool asynchronously."""
        raise NotImplementedError()
