from adme_pred import ADME
from langchain.tools import BaseTool
from rdkit import Chem


class CalculateBBBPermeant(BaseTool):
    name = "CalculateBBBPermeant"
    description = "calculates the Blood Brain Barrier Permeability of the compound"

    def _run(self, compound_smiles: str) -> str:
        """
        Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master
            From the adme-pred-py github:
        Daina (2016) A BOILED-Egg To Predict Gastrointestinal Absorption and
            Brain Penetration of Small Molecules

            This multivariate model uses log P and Polar Surface Area to determine
            druglikeness. This function implements their Blood Brain Barrier
            (BBB) model, which is the "yolk" of the BOILED-Egg.
        """

        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        if mol.boiled_egg_bbb():
            return "Yes"
        else:
            return "No"

    async def _arun(self, compound_smiles: str) -> str:
        """Use the calculate_bbb_permeant tool asynchronously."""
        raise NotImplementedError()


# example use
# sample = calculate_bbb_permeant("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
# print(sample)
