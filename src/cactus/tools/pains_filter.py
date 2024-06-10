"""Tool for calculating if a molecule passes the PAINS Filter."""

from adme_pred import ADME
from langchain.tools import BaseTool
from rdkit import Chem

DESC = """
Used when you need to calculate whether a molecule triggers the Pains Filter.
"""


class PainsFilter(BaseTool):
    """Calculates the PAINS Filter."""

    name = "PainsFilter"
    description = DESC

    def _run(self, compound_smiles: str) -> bool:
        """Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.

        From the adme-pred-py github:
        Baell and Holloway (2010) New Substructure Filters for Removal of Pan
        Assay Interference Compounds (PAINS) from Screening Libraries and for
        Their Exclusion in Bioassays

        This filter finds promiscuous compounds that are likely to show activity
        regardless of the target.

        Returns
        -------
            Boolean of whether the molecule triggers the PAINS filter.
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        return mol.pains()

    async def _arun(self, compound_smiles: str) -> bool:
        """Use the brenk_filter tool asynchronously."""
        raise NotImplementedError()
