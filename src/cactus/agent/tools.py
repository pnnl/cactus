"""Script for loading tools for the Cactus Agent."""
from langchain import agents

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
    CasToSMILES,
    ChemblidToSMILES,
    CidToSMILES,
    InchikeyToSMILES,
    MolecularFormulaToSMILES,
    NameToSMILES,
    PainsFilter,
    ZincIDToSMILES,
)


def make_tools():
    """Method for aggregating and generating a list of tools for the LLM Agent"""
    all_tools = agents.load_tools(
        [
            "ddg-search",
            "wikipedia",
        ]
    )

    all_tools += [
        InchikeyToSMILES(),
        NameToSMILES(),
        # CasToSMILES(),
        # ChemblidToSMILES(),
        CidToSMILES(),
        # MolecularFormulaToSMILES(),
        # ZincIDToSMILES(),
        CalculateMolWt(),
        CalculateQED(),
        BrenkFilter(),
        CalculateTPSA(),
        # CalculateBBBPermeant(),
        # CalculateDruglikeness(),
        # CalculateGIAbsorption(),
        CalculateLogP(),
        # PainsFilter(),
        CalculateSA(),
    ]

    return all_tools
