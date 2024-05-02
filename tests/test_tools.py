import pytest

from cactus.tools import (
    BrenkFilter,
    CalculateBBBPermeant,
    CalculateDruglikeness,
    CalculateGIAbsorption,
    CalculateLogP,
    CalculateMolWt,
    CalculateQED,
    CalculateSA,
    CalculateTPSA,
    PainsFilter,
)


# Brenk Filter Test
def test_brenk_filter_success():
    """Tests BrenkFilter tool with a valid SMILES string."""
    compound_smiles = "CC(=O)N"
    tool = BrenkFilter()
    brenk_pass_fail = tool.run(compound_smiles)
    assert brenk_pass_fail == False  # Brenk should either be True or False


def test_calculate_qed_success():
    """Tests CalculateQED tool with a valid SMILES string."""
    compound_smiles = "CC(=O)N"
    tool = CalculateQED()
    qed_score = tool.run(compound_smiles)
    assert 0 <= qed_score <= 1  # QED score should be between 0 and 1
    # You can optionally compare the score with an expected value from RDKit


def test_calculate_qed_invalid_smiles():
    """Tests CalculateQED tool with an invalid SMILES string."""
    invalid_smiles = "invalid_smiles"
    tool = CalculateQED()
    with pytest.raises(ValueError):
        tool.run(invalid_smiles)
