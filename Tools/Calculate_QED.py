from rdkit.Chem import MolFromSmiles, Descriptors

def calculate_QED(compound: str) -> float:
    """
    Compute quantum estimate of druglikeness (QED) of the given molecule.
    
    Parameters:
    compound: Compound in SMILES format

    Returns:
    float: The QED (0 to 1)
    """
    return Descriptors.qed(MolFromSmiles(compound))