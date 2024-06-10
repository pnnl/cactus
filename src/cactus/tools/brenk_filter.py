"""Tool for calculating if a compound passes the Brenk Filter."""

from adme_pred import ADME
from langchain.tools import BaseTool
from rdkit import Chem

DESC = """
Used when you need to calculate whether a molecule triggers the Brenk Filter.
"""


class BrenkFilter(BaseTool):
    """Tool to check if the molecule passes the Brenk Filter."""

    name = "BrenkFilter"
    description = DESC

    def _run(self, compound_smiles: str) -> bool:
        """Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.

            From the adme-pred-py documentation:
        Brenk (2008) Lessons Learnt from Assembling Screening Libraries for
            Drug Discovery for Neglected Diseases

            Brenk's Structural Alert filter finds fragments "putatively toxic,
            chemically reactive, metabolically unstable or to bear properties
            responsible for poor pharmacokinetics."

        Returns
        -------
                Boolean of whether the molecule triggers the Brenk filter.
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        return mol.brenk()

    async def _arun(self, compound_smiles: str) -> bool:
        """Use the brenk_filter tool asynchronously."""
        raise NotImplementedError()
