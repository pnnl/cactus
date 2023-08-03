from rdkit.Chem import MolFromSmiles, Descriptors

def calculate_TPSA(compound: str) -> float:
    """
    Compute the topological surface area (TPSA) of the given molecule.
    
    Parameters:
    compound: Compound in SMILES format

    Returns:
    float: The TPSA in angstroms^2
    """
    return Descriptors.TPSA(MolFromSmiles(compound))