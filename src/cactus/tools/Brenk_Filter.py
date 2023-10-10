from langchain.tools import BaseTool
from rdkit import Chem
from adme_pred import ADME


class brenk_filter(BaseTool):
    name = "brenk_filter"
    description = "returns a boolean value of whether the molecule triggers the Brenk Filter"

    def _run(self, compound_smiles: str) -> bool:
        """
        Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master

            From the adme-pred-py documentation:
        Brenk (2008) Lessons Learnt from Assembling Screening Libraries for
            Drug Discovery for Neglected Diseases

            Brenk's Structural Alert filter finds fragments "putatively toxic,
            chemically reactive, metabolically unstable or to bear properties
            responsible for poor pharmacokinetics."

            Returns:
                Boolean of whether the molecule triggers the Brenk filter.
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        return mol.brenk()

    async def _arun(self, compound_smiles: str) -> bool:
        """Use the brenk_filter tool asynchronously."""
        raise NotImplementedError()


# example use
# sample = brenk_filter("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
# print(sample)
