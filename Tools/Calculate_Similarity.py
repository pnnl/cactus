from rdkit.Chem import MolFromSmiles, AllChem, MACCSkeys, RDKFingerprint
from scipy.spatial import distance

## This should be the default
def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Jaccard-Tanimoto) similarity value between exactly 2 molecules (Morgan).
    
    These papers expand on why this should be default:
    Jaccard - https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0069-3#Sec12
    Morgan  - https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0148-0#Sec14

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Jaccard similarity between the two compound's ECFP6 (Morgan) fingerprints.
    """
    try: compound_1 = AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(compound_1),3,nBits=1024)
    except: return "Compound 1 is invalid."
    try: compound_2 = AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(compound_2),3,nBits=1024)
    except: return "Compound 2 is invalid."

    return 1 - distance.jaccard(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

## Other Variations:

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (Morgan).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Dice similarity between the two compound's ECFP6 (Morgan) fingerprints.
    """
    try: compound_1 = AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(compound_1),3,nBits=1024)
    except: return "Compound 1 is invalid."
    try: compound_2 = AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(compound_2),3,nBits=1024)
    except: return "Compound 2 is invalid."

    return 1 - distance.dice(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Cosine) similarity value between exactly 2 molecules (Morgan).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Cosine similarity between the two compound's ECFP6 (Morgan) fingerprints.
    """
    try: compound_1 = AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(compound_1),3,nBits=1024)
    except: return "Compound 1 is invalid."
    try: compound_2 = AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(compound_2),3,nBits=1024)
    except: return "Compound 2 is invalid."

    return 1 - distance.cosine(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Cosine) similarity value between exactly 2 molecules (MACCS).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Cosine similarity between the two compound's (MACCS) fingerprints.
    """
    try: compound_1 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_1))
    except: return "Compound 1 is invalid."
    try: compound_2 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_2))
    except: return "Compound 2 is invalid."

    return 1 - distance.cosine(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (MACCS).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Dice similarity between the two compound's (MACCS) fingerprints.
    """
    try: compound_1 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_1))
    except: return "Compound 1 is invalid."
    try: compound_2 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_2))
    except: return "Compound 2 is invalid."

    return 1 - distance.dice(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Jaccard) similarity value between exactly 2 molecules (MACCS).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Jaccard similarity between the two compound's (MACCS) fingerprints.
    """
    try: compound_1 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_1))
    except: return "Compound 1 is invalid."
    try: compound_2 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_2))
    except: return "Compound 2 is invalid."

    return 1 - distance.jaccard(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Jaccard) similarity value between exactly 2 molecules (Daylight).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Jaccard similarity between the two compound's (RDKit/Daylight) fingerprints.
    """
    try: compound_1 = RDKFingerprint(MolFromSmiles(compound_1))
    except: return "Compound 1 is invalid."
    try: compound_2 = RDKFingerprint(MolFromSmiles(compound_2))
    except: return "Compound 2 is invalid."

    return 1 - distance.jaccard(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (Daylight).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Dice similarity between the two compound's (RDKit/Daylight) fingerprints.
    """
    try: compound_1 = RDKFingerprint(MolFromSmiles(compound_1))
    except: return "Compound 1 is invalid."
    try: compound_2 = RDKFingerprint(MolFromSmiles(compound_2))
    except: return "Compound 2 is invalid."

    return 1 - distance.dice(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity

def calculate_similarity(compound_1: str, compound_2: str) -> float:
    """
    Calculate (Cosine) similarity value between exactly 2 molecules (Daylight).

    Parameters:
    compound_1: Compound in SMILES format
    compound_2: Compound in SMILES format

    Returns:
    float: Cosine similarity between the two compound's (RDKit/Daylight) fingerprints.
    """
    try: compound_1 = RDKFingerprint(MolFromSmiles(compound_1))
    except: return "Compound 1 is invalid."
    try: compound_2 = RDKFingerprint(MolFromSmiles(compound_2))
    except: return "Compound 2 is invalid."

    return 1 - distance.cosine(compound_1,compound_2) # since scipy calculates distance, 1-d = similarity