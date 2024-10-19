"""Tool to calculate if a compound passes the Lipinski rule of five."""

from adme_py import ADME
from langchain.tools import BaseTool


class CalculateDruglikeness(BaseTool):
    """Tool to calculate if a compound passes the Lipinski Rule of Five."""

    name: str = "calculate_druglikeness"
    description: str = (
        "calculates the druglikeness of the compound with regards to Lipinski's rule of 5"
    )

    def _run(self, compound_smiles: str) -> str:
        """Calculate the lipinski druglikeness of the compound.

        Parameters
        ----------
        compound_smiles : str
            The input smiles string.

        Returns
        -------
        Union[str, dict[str, str]]
            - "Pass" if the molecule adheres to all rules.
            - A dictionary with the violated rules as keys and descriptive messages as values.

        Notes
        -----
        The Lipinksi Rule of 5 is described in:
        C.A Lipinksi, et al. Adv. Drug Delivery Rev. 2001, 46, 1-3, 3-26
        https://doi.org/10.1016/S0169-409X(00)00129-0
        """
        summary = ADME(compound_smiles).calculate()
        return summary["druglikeness"]["lipinski"]

    async def _arun(self, compound_smiles: str) -> str:
        """Use the calculate_druglikeness tool asynchronously."""
        raise NotImplementedError()
