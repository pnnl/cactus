from rdkit.Chem import MolFromSmiles, AllChem, MACCSkeys, RDKFingerprint
from scipy.spatial import distance

## Helper function to fix some SMILES inconsistencies instead of always spitting it out:

from rdkit.Chem import SanitizeMol, SanitizeFlags
def partially_sanitize(smi):
    """
    From the RDKit Documentation: https://www.rdkit.org/docs/Cookbook.html#explicit-valence-error-partial-sanitization

    Mostly just makes sure the charge of any hypervalent atoms is explicit for RDKit. (In my experience, very common in Nitrogen, Oxygen, and Sulfur)
    I usually need this when working with Fingerprints (they are very strict)

    I'm noticing some weird behavior with ignoring SMILES' implicit hydrogens but I need to test more. 
    For now, just using MolFromSmiles() and no OPENBABEL or manually remove implicit hypervalence is fine
    e.g. n implicit isonitriles, nitrates, or quaternary amines (CN#C -> C[N+]#[C-])

    I just think this is necessary for the final tool to ingest *any* smiles or slightly wrong ones (see: LLM)
    ESPECIALLY if we're also using RDKit or ML to generate new molecule strings
    """
    if '\t' in smi: smi = smi.split('\t')[0] # weird OPENBABEL syntax
    mol = MolFromSmiles(smi,sanitize=False) # if True, this is where RDKit throws error. False=skip valence calc
    mol.UpdatePropertyCache(strict=False)  # recalculate the valences, be lenient
    # then, redo the molecule with charges, aromaticity, etc explicitly stated
    SanitizeMol(mol,SanitizeFlags.SANITIZE_FINDRADICALS|SanitizeFlags.SANITIZE_KEKULIZE|SanitizeFlags.SANITIZE_SETAROMATICITY|SanitizeFlags.SANITIZE_SETCONJUGATION|SanitizeFlags.SANITIZE_SETHYBRIDIZATION|SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
    return mol

## This should be the default:

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