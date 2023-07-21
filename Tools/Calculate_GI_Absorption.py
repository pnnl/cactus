from rdkit import Chem
from adme_pred import ADME

def calculate_gi_absorption(compound_smiles: str) -> str:
	"""
	Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master
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

#example use
sample = calculate_gi_absorption("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
print(sample)

