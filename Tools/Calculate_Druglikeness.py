from rdkit import Chem
from adme_pred import ADME

def calculate_druglikeness(compound_smiles: str) -> str:
	"""
	Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master

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

#example use
sample = calculate_druglikeness("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
print(sample)
