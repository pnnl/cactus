from rdkit.Chem import MolFromSmiles, Descriptors

def calculate_QED(compound: str) -> float:
    """
    Compute Quantitative Estimate of Druglikeness (QED) of the given molecule. Bickerton et al 2012.
    
    Parameters:
    compound: Compound in SMILES format

    Returns:
    float: The QED from 0 (druglike) to 1 (not)
    """
    return Descriptors.qed(MolFromSmiles(compound))