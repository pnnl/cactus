"""Method script for generating some aspects of the chem benchmark"""
from typing import List

import itertools
import pandas as pd


def generate_permutations(template: str, attribute: List, compounds: List):
    """Generate csv file of questions from list of compounds and descriptors"""
    questions = pd.DataFrame(columns=["question"])
    new_rows = []
    all_permutations = list(itertools.product(attribute, compounds))

    for permutation in all_permutations:
        new_row = {"question": template.format(attribute=permutation[0], compound=permutation[1])}
        new_rows.append(new_row)

    questions = pd.concat([questions, pd.DataFrame(new_rows)], ignore_index=True)

    questions.to_csv("QuestionsChem.csv")


if __name__ == "__main__":
    TEMPLATE = "What is the {attribute} of {compound}?"
    pubchem_molecules = [
        "Acetaminophen",
        "Aspirin",
        "Caffeine",
        "Ethanol",
        "Glucose",
        "Hydrocodone",
        "Methanol",
        "Morphine",
        "Nicotine",
        "Oxygen",
        "Penicillin",
        "Water",
        "Vitamin C",
        "Acetic acid",
        "Ammonia",
        "Carbon dioxide",
        "Ethyl acetate",
        "Glucose 6-phosphate",
        "Methane",
        "Sulfur dioxide",
        "Sucrose",
        "Benzaldehyde",
        "CC(=O)NC1=CC=C(C=C1)O",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CCO",
        "C(C1C(C(C(C(O1)O)O)O)O)O",
        "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4",
        "CO",
        "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        "CN1CCCC1C2=CN=CC=C2",
        "O=O",
        "CC1(C(N2C(S1=O)C(C2=O)NC(=O)C3=CC=CC=C3)C(=O)OC(C4=CC=CC=C4)C5=CC=CC=C5)C",
        "O",
        "C(C(C1C(=C(C(=O)O1)O)O)O)O",
        "CC(=O)O",
        "N",
        "C(=O)=O",
        "CCOC(=O)C",
        "C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O",
        "C",
        "O=S=O",
        "C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
        "C1=CC=C(C=C1)C=O",
    ]

    descriptors = [
        "Molecular Weight",
        "QED",
        "Synthetic Accesibility",
        "Total Polar Surface Area",
        "LogP",
    ]

    generate_permutations(TEMPLATE, descriptors, pubchem_molecules)
