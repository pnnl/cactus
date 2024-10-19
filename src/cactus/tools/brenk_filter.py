"""Tool for calculating if a compound passes the Brenk Filter."""

from adme_py import ADME
from langchain.tools import BaseTool

DESC = """
Used when you need to calculate whether a molecule triggers the Brenk Filter.
"""


class BrenkFilter(BaseTool):
    """Tool to check if the molecule passes the Brenk Filter."""

    name: str = "BrenkFilter"
    description: str = DESC

    def _run(self, compound_smiles: str) -> bool:
        """Check if a molecule triggers the Brenk filter.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The input RDKit molecule object.

        Returns
        -------
        bool
            True if the molecule triggers the Brenk filter, False otherwise.

        Notes
        -----
        The Brenk filter is described in:
        R. Brenk, et al., ChemMedChem, 2008, 3: 435-444.
        https://doi.org/10.1002/cmdc.200700139
        """
        summary = ADME(compound_smiles).calculate()

        result = summary["medicinal"]["brenk"]

        if result:
            return "Does not pass, triggers the Brenk Filter."
        else:
            return "Passes, does not trigger the Brenk Filter."

    async def _arun(self, compound_smiles: str) -> bool:
        """Use the brenk_filter tool asynchronously."""
        raise NotImplementedError()
