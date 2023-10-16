import os
from langchain import agents
from langchain.base_language import BaseLanguageModel
from cactus.tools import (
    cas_to_SMILES,
    chemblid_to_SMILES,
    cid_to_SMILES,
    molecular_formula_to_SMILES,
    inchikey_to_SMILES,
    name_to_SMILES,
    zinc_id_to_SMILES,
    calculate_MolWt,
    calculate_QED,
    calculate_TPSA,
    calculate_logp,
    calculate_SA,
    calculate_bbb_permeant,
    pains_filter,
    calculate_gi_absorption,
    brenk_filter,
    calculate_druglikeness,
)


def make_tools():
    """Method for aggregating and generating a list of tools for the LLM Agent"""
    all_tools = agents.load_tools(
        [
            "python_repl",
            "ddg-search",
            "wikipedia",
        ]
    )

    all_tools += [
        inchikey_to_SMILES(),
        name_to_SMILES(),
        cas_to_SMILES(),
        chemblid_to_SMILES(),
        cid_to_SMILES(),
        molecular_formula_to_SMILES(),
        zinc_id_to_SMILES(),
        calculate_MolWt(),
        calculate_QED(),
        brenk_filter(),
        calculate_TPSA(),
        calculate_bbb_permeant(),
        calculate_druglikeness(),
        calculate_gi_absorption(),
        calculate_logp(),
        pains_filter(),
        calculate_SA(),
    ]

    return all_tools
