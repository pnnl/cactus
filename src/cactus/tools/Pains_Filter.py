from rdkit import Chem
from adme_pred import ADME

<<<<<<< HEAD

def pains_filter(self, compound_smiles: str) -> bool:
    """
    Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master
    From the adme-pred-py github:
    Baell and Holloway (2010) New Substructure Filters for Removal of Pan
    Assay Interference Compounds (PAINS) from Screening Libraries and for
    Their Exclusion in Bioassays

    This filter finds promiscuous compounds that are likely to show activity
    regardless of the target.

    Returns:
        Boolean of whether the molecule triggers the PAINS filter.
    """
    mol = Chem.MolFromSmiles(compound_smiles)
    mol = ADME(mol)
    return mol.pains()


# example use
=======
def pains_filter(compound_smiles: str) -> bool:
	"""
	Uses the adme-pred-py implementation: https://github.com/ikmckenz/adme-pred-py/tree/master
        From the adme-pred-py github:
	Baell and Holloway (2010) New Substructure Filters for Removal of Pan
        Assay Interference Compounds (PAINS) from Screening Libraries and for
        Their Exclusion in Bioassays

        This filter finds promiscuous compounds that are likely to show activity
        regardless of the target.

        Returns:
            Boolean of whether the molecule triggers the PAINS filter.
	"""
	mol = Chem.MolFromSmiles(compound_smiles)
	mol = ADME(mol)
	return mol.pains()

#example use
>>>>>>> 95cc70a (refactoring to src/cactus)
sample = pains_filter("CC(Cc1ccc(cc1)C(C(=O)O)C)C")
print(sample)
