import os
import sys
from rdkit.Chem import MolFromSmiles, RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

def calculate_SA(compound: str) -> float:
    """
    Compute Synthetic Accessibility (SA) of the given molecule. Ertl & Schuffenhauer 2009.
    
    Parameters:
    compound: Compound in SMILES format

    Returns:
    float: The SA between 1 (easy) and 10 (hard)
    """
    return sascorer.calculateScore(MolFromSmiles(compound))