from rdkit.Chem import MolFromSmiles, Descriptors

def calculate_MolWt(compound: str) -> float:
    """
    Compute the exact molecular weight of the given molecule.
    
    Parameters:
    compound: Compound in SMILES format

    Returns:
    float: The exact molecular weight (0 to ∞)
    """
    return Descriptors.ExactMolWt(MolFromSmiles(compound))