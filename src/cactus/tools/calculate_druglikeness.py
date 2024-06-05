"""Tool to calculate if a compound passes the Lipinski rule of five."""

from adme_pred import ADME
from langchain.tools import BaseTool
from rdkit import Chem


class CalculateDruglikeness(BaseTool):
    """Tool to calculate if a compound passes the Lipinski Rule of Five."""

    name = "calculate_druglikeness"
    description = "calculates the druglikeness of the compound with regards to Lipinski's rule of 5"

    def _run(self, compound_smiles: str) -> str:
        """Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.

        From the adme-pred-py github:
        Lipinski (2001) Experimental and computational approaches to estimate
            solubility and permeability in drug discovery and development settings

            https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five

            Lipinski's rule of 5 is one of the most important druglikenss filters,
            against which all others are judged. The rules of the filter are no
            more than 5 hydrogen bond donors, no more than 10 hydrogen bond
            acceptors, a molecular mass of less than 500 daltons, and a logP
            that does not exceed 5.
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        return mol.druglikeness_lipinski(verbose=True)

    async def _arun(self, compound_smiles: str) -> str:
        """Use the calculate_druglikeness tool asynchronously."""
        raise NotImplementedError()
