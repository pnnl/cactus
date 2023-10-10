from rdkit import Chem
from adme_pred import ADME
from rdkit.Chem import Descriptors

<<<<<<< HEAD

# *******************Brenk***************************
def brenk_filter(compound_smiles: str) -> bool:
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


# *****************Lipinski druglikeness check************
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


# **************Blood Brain Barrier Permeant Calculation **********
def calculate_bbb_permeant(compound_smiles: str) -> str:
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


# *************** GI Absorption Calculation ******************
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


# *************** LogP Calculation ***********************
def calculate_logp(compound_smiles: str) -> float:
    """
    RdKit descriptor module to calculate LogP from the atom based Crippen approach
    https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
    https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
    """
    mol = Chem.MolFromSmiles(compound_smiles)
    return Descriptors.MolLogP(mol)


# ***************PAINS Check****************************
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


# example use
test_smiles = "CC(Cc1ccc(cc1)C(C(=O)O)C)C"
=======
"""
*******************Brenk***************************
"""
def brenk_filter(compound_smiles: str) -> bool:
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

"""
*****************Lipinski druglikeness check************
"""
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

"""
**************Blood Brain Barrier Permeant Calculation **********
"""
def calculate_bbb_permeant(compound_smiles: str) -> str:
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

"""
*************** GI Absorption Calculation ******************
"""
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

"""
*************** LogP Calculation ***********************
"""
def calculate_logp(compound_smiles: str) -> float:
	"""
        RdKit descriptor module to calculate LogP from the atom based Crippen approach
        https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
        https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html
	"""
	mol = Chem.MolFromSmiles(compound_smiles)
	return Descriptors.MolLogP(mol)

"""
***************PAINS Check****************************
"""
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
test_smiles ="CC(Cc1ccc(cc1)C(C(=O)O)C)C"
>>>>>>> 95cc70a (refactoring to src/cactus)

sample = brenk_filter(test_smiles)
print("Brenk Filter: ", sample)

sample = calculate_druglikeness(test_smiles)
print("Lipinski: ", sample)

sample = calculate_bbb_permeant(test_smiles)
print("BBB: ", sample)

sample = calculate_gi_absorption(test_smiles)
print("GI: ", sample)

sample = calculate_logp(test_smiles)
print("LogP: ", sample)

sample = pains_filter(test_smiles)
print("PAINS Filter: ", sample)
<<<<<<< HEAD
=======

>>>>>>> 95cc70a (refactoring to src/cactus)
