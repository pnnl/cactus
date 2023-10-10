from langchain.tools import BaseTool
from rdkit import Chem
from rdkit.Chem import Descriptors


class calculate_logp(BaseTool):
    name = "calculate_logp"
    description = "Returns a Boolean of whether the molecule triggers the PAINS filter"

    def _run(self, compound_smiles: str) -> float:
        """
        RdKit descriptor module to calculate LogP from the atom based Crippen approach
        https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
        https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
        """
        mol = Chem.MolFromSmiles(compound_smiles)
        return Descriptors.MolLogP(mol)

    async def _arun(self, compound_smiles: str) -> float:
        """Use the calculate_logp tool asynchronously."""
        raise NotImplementedError()

    # example use


# sample = calculate_logp("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
# print(sample)
