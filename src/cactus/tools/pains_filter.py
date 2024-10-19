"""Tool for calculating if a molecule passes the PAINS Filter."""

from adme_py import ADME
from langchain.tools import BaseTool

DESC = """
Used when you need to calculate whether a molecule triggers the PAINS Filter.
"""


class PainsFilter(BaseTool):
    """Calculates the PAINS Filter."""

    name: str = "PainsFilter"
    description: str = DESC

    def _run(self, compound_smiles: str) -> bool:
        """Check if a molecule triggers the PAINS (Pan Assay Interference Compounds) filter.

        Parameters
        ----------
        compound_smiles : str
            The input smiles for the molecule.

        Returns
        -------
        bool
            True if the molecule triggers the PAINS filter, False otherwise.

        Notes
        -----
        The PAINS filter is described in:
        J.B. Baell and G.A. Holloway, J. Med. Chem. 2010, 53, 7, 2719–2740
        https://doi.org/10.1021/jm901137j
        """
        summary = ADME(compound_smiles).calculate()

        result = summary["medicinal"]["pains"]

        if result:
            return "Does not pass, triggers the PAINS Filter."
        else:
            return "Passes, does not trigger the PAINS Filter."

    async def _arun(self, compound_smiles: str) -> bool:
        """Use the pains_filter tool asynchronously."""
        raise NotImplementedError()
