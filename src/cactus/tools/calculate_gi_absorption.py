"""Tool to calculate the GI Absorption of the compound."""

from adme_pred import ADME
from langchain.tools import BaseTool
from rdkit import Chem


class CalculateGIAbsorption(BaseTool):
    """Calculate the GI Absorption of a compound."""

    name = "calculate_gi_absorption"
    description = "returns whether the gastrointestinal absorption is high or low"

    def _run(self, compound_smiles: str) -> str:
        """Use the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master.

        From the adme-pred-py github:
            Daina (2016) A BOILED-Egg To Predict Gastrointestinal Absorption and
            Brain Penetration of Small Molecules

            This multivariate model uses log P and Polar Surface Area to determine
            druglikeness. This function implements their Human Intestinal Absorption
            (HIA) model, which is the "white" of the BOILED-Egg.
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        mol = ADME(mol)
        if mol.boiled_egg_hia():
            return "High"
        else:
            return "Low"

    async def _arun(self, compound_smiles: str) -> str:
        """Use the calculate_gi_absorption tool asynchronously."""
        raise NotImplementedError()
