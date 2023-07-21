from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_logp(compound_smiles: str) -> float:
	"""
	RdKit descriptor module to calculate LogP from the atom based Crippen approach
	https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
	https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
	"""
	mol = Chem.MolFromSmiles(compound_smiles)
	return Descriptors.MolLogP(mol)

#example use
sample = calculate_logp("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
print(sample)
