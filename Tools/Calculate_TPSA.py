from rdkit.Chem import MolFromSmiles, Descriptors

def calculate_TPSA(compound: str) -> float:
    """
    Compute the topological surface area (TPSA) of the given molecule.
    
    Parameters:
    compound: Compound in SMILES format

    Returns:
    float: The TPSA
    """
    return Descriptors.TPSA(MolFromSmiles(compound))