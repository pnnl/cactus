# ruff: noqa
"""Tools to Calculate Molecular Similarity."""

from langchain.tools import BaseTool
from rdkit.Chem import (
    AllChem,
    MACCSkeys,
    MolFromSmiles,
    RDKFingerprint,
    SanitizeFlags,
    SanitizeMol,
    DataStrucs,
)
from scipy.spatial import distance
import numpy as np


class PartiallySanitize(BaseTool):
    """Run initial sanitize over compound."""

    name = "partially_sanitize"
    description = "To check for explicit valence error by doing a partial sanitization"

    def _run(self, smi: str) -> str:
        """Ensure charge of hypervalent atoms are explicitly for RDKit.

        From the RDKit Documentation: https://www.rdkit.org/docs/Cookbook.html#explicit-valence-error-partial-sanitization
        """
        if "\t" in smi:
            smi = smi.split("\t")[0]  # weird OPENBABEL syntax
        mol = MolFromSmiles(
            smi, sanitize=False
        )  # if True, this is where RDKit throws error. False=skip valence calc
        mol.UpdatePropertyCache(strict=False)  # recalculate the valences, be lenient
        # then, redo the molecule with charges, aromaticity, etc explicitly stated
        SanitizeMol(
            mol,
            SanitizeFlags.SANITIZE_FINDRADICALS
            | SanitizeFlags.SANITIZE_KEKULIZE
            | SanitizeFlags.SANITIZE_SETAROMATICITY
            | SanitizeFlags.SANITIZE_SETCONJUGATION
            | SanitizeFlags.SANITIZE_SETHYBRIDIZATION
            | SanitizeFlags.SANITIZE_SYMMRINGS,
            catchErrors=True,
        )
        return mol

    async def _arun(self, smi: str) -> str:
        """Use the partially_sanitize tool asynchronously."""
        raise NotImplementedError()


# This should be the default:
class CalculateJaccardSimilarity(BaseTool):
    """Calculate the Jaccard similarity between 2 compounds."""

    name = "calculate_similarity"
    description = (
        "Calculate (Jaccard-Tanimoto) similarity value between exactly 2 molecules (Morgan)."
    )

    def _run(self, compound_1: str, compound_2: str) -> float:
        """Calculate (Jaccard-Tanimoto) similarity value between exactly 2 molecules (Morgan).

        These papers expand on why this should be default:
        Jaccard - https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0069-3#Sec12
        Morgan  - https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0148-0#Sec14

        Parameters
        ----------
             compound_1: Compound in SMILES format
             compound_2: Compound in SMILES format

        Returns
        -------
             float: Jaccard similarity between the two compound's ECFP6 (Morgan) fingerprints.
        """
        try:
            mol1 = MolFromSmiles(compound_1)

            if mol1 is None:
                raise ValueError("Compound 1 is invalid.")

            compound_1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=1024)
        except ValueError:
            return np.nan

        try:
            mol2 = MolFromSmiles(compound_2)

            if mol2 is None:
                raise ValueError("Compound 2 is invalid.")

            compound_2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=1024)

        except ValueError:
            return np.nan

        return 1 - distance.jaccard(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    # Other Variations:
    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class CalculateDiceSimilarity(BaseTool):
    """Calculate the Sorensen-Dice similarity."""

    name = "calculate_similarity"
    description = "Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (Morgan)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (Morgan).

        Parameters
        ----------
             compound_1: Compound in SMILES format
             compound_2: Compound in SMILES format

        Returns
        -------
             float: Dice similarity between the two compound's ECFP6 (Morgan) fingerprints.
        """
        try:
            mol1 = MolFromSmiles(compound_1)

            if mol1 is None:
                raise ValueError("Compound 1 is invalid.")

            compound_1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=1024)
        except ValueError:
            return np.nan

        try:
            mol2 = MolFromSmiles(compound_2)

            if mol2 is None:
                raise ValueError("Compound 2 is invalid.")

            compound_2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=1024)

        except ValueError:
            return np.nan

        return 1 - distance.dice(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class CalculateCosineSimilarity(BaseTool):
    """Calculate the Cosine Similarity between 2 compounds."""

    name = "calculate_similarity"
    description = "Calculate (Cosine) similarity value between exactly 2 molecules (Morgan)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """Calculate (Cosine) similarity value between exactly 2 molecules (Morgan).

        Parameters
        ----------
             compound_1: Compound in SMILES format
             compound_2: Compound in SMILES format

        Returns
        -------
             float: Cosine similarity between the two compound's ECFP6 (Morgan) fingerprints.
        """
        try:
            mol1 = MolFromSmiles(compound_1)

            if mol1 is None:
                raise ValueError("Compound 1 is invalid.")

            compound_1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=1024)
        except ValueError:
            return np.nan

        try:
            mol2 = MolFromSmiles(compound_2)

            if mol2 is None:
                raise ValueError("Compound 2 is invalid.")

            compound_2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=1024)

        except ValueError:
            return np.nan

        return 1 - distance.cosine(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class CalculateCosineSimilarityMaccs(BaseTool):
    """Calculate the cosine similarity between 2 compounds using MACCS."""

    name = "calculate_similarity"
    description = "Calculate (Cosine) similarity value between exactly 2 molecules (MACCS)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """Calculate (Cosine) similarity value between exactly 2 molecules (MACCS).

        Parameters
        ----------
             compound_1: Compound in SMILES format
             compound_2: Compound in SMILES format

        Returns
        -------
             float: Cosine similarity between the two compound's (MACCS) fingerprints.
        """
        try:
            compound_1 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_1))
        except:
            return "Compound 1 is invalid."
        try:
            compound_2 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_2))
        except:
            return "Compound 2 is invalid."

        return 1 - distance.cosine(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class calculate_similarity(BaseTool):
    name = "calculate_similarity"
    description = "Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (MACCS)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """
        Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (MACCS).

        Parameters
        ----------
        compound_1: Compound in SMILES format
        compound_2: Compound in SMILES format

        Returns
        -------
        float: Dice similarity between the two compound's (MACCS) fingerprints.
        """
        try:
            compound_1 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_1))
        except:
            return "Compound 1 is invalid."
        try:
            compound_2 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_2))
        except:
            return "Compound 2 is invalid."

        return 1 - distance.dice(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class calculate_similarity(BaseTool):
    name = "calculate_similarity"
    description = " Calculate (Jaccard) similarity value between exactly 2 molecules (Daylight)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """
        Calculate (Jaccard) similarity value between exactly 2 molecules (MACCS).

        Parameters
        ----------
        compound_1: Compound in SMILES format
        compound_2: Compound in SMILES format

        Returns
        -------
        float: Jaccard similarity between the two compound's (MACCS) fingerprints.
        """
        try:
            compound_1 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_1))
        except:
            return "Compound 1 is invalid."
        try:
            compound_2 = MACCSkeys.GenMACCSKeys(MolFromSmiles(compound_2))
        except:
            return "Compound 2 is invalid."

        return 1 - distance.jaccard(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class calculate_similarity(BaseTool):
    name = "calculate_similarity"
    description = " Calculate (Jaccard) similarity value between exactly 2 molecules (Daylight)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """
        Calculate (Jaccard) similarity value between exactly 2 molecules (Daylight).

        Parameters
        ----------
        compound_1: Compound in SMILES format
        compound_2: Compound in SMILES format

        Returns
        -------
        float: Jaccard similarity between the two compound's (RDKit/Daylight) fingerprints.
        """
        try:
            compound_1 = RDKFingerprint(MolFromSmiles(compound_1))
        except:
            return "Compound 1 is invalid."
        try:
            compound_2 = RDKFingerprint(MolFromSmiles(compound_2))
        except:
            return "Compound 2 is invalid."

        return 1 - distance.jaccard(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class calculate_similarity(BaseTool):
    name = "calculate_similarity"
    description = (
        "Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (Daylight)"
    )

    def _run(self, compound_1: str, compound_2: str) -> float:
        """
        Calculate (Sørensen-Dice) similarity value between exactly 2 molecules (Daylight).

        Parameters
        ----------
        compound_1: Compound in SMILES format
        compound_2: Compound in SMILES format

        Returns
        -------
        float: Dice similarity between the two compound's (RDKit/Daylight) fingerprints.
        """
        try:
            compound_1 = RDKFingerprint(MolFromSmiles(compound_1))
        except:
            return "Compound 1 is invalid."
        try:
            compound_2 = RDKFingerprint(MolFromSmiles(compound_2))
        except:
            return "Compound 2 is invalid."

        return 1 - distance.dice(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()


class CalculateSimilarity(BaseTool):
    name = "calculate_similarity"
    description = "Calculate (Cosine) similarity value between exactly 2 molecules (Daylight)."

    def _run(self, compound_1: str, compound_2: str) -> float:
        """
        Calculate (Cosine) similarity value between exactly 2 molecules (Daylight).

        Parameters
        ----------
        compound_1: Compound in SMILES format
        compound_2: Compound in SMILES format

        Returns
        -------
        float: Cosine similarity between the two compound's (RDKit/Daylight) fingerprints.
        """
        try:
            compound_1 = RDKFingerprint(MolFromSmiles(compound_1))
        except:
            return "Compound 1 is invalid."
        try:
            compound_2 = RDKFingerprint(MolFromSmiles(compound_2))
        except:
            return "Compound 2 is invalid."

        return 1 - distance.cosine(
            compound_1, compound_2
        )  # since scipy calculates distance, 1-d = similarity

    async def _arun(self, compound_1: str, compound_2: str) -> float:
        """Use the calculate_similarity tool asynchronously."""
        raise NotImplementedError()
