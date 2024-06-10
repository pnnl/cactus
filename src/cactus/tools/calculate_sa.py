"""Tool to calculate the SA of a compound."""

import os
import sys

from rdkit.Chem import MolFromSmiles, RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer
from langchain.tools import BaseTool


class CalculateSA(BaseTool):
    """Calculate the SA of the compound."""

    name = "CalculateSA"
    description = "Used to compute the synthetic accessibility (SA) of the given molecule."

    def _run(self, compound: str) -> float:
        """Compute Synthetic Accessibility (SA) of the given SMILES string. Ertl & Schuffenhauer 2009.

        Parameters
        ----------
        compound: Compound in SMILES format

        Returns
        -------
        float: The SA between 1 (easy) and 10 (hard)
        """
        return sascorer.calculateScore(MolFromSmiles(compound))

    async def _arun(self, compound: str) -> float:
        """Use the convert_to_SMILES tool asynchronously."""
        raise NotImplementedError()
