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
        "Insulin",
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
    ]

    descriptors = [
        "Molecular Weight",
        "QED",
        "Synthetic Accesibility",
        "Total Polar Surface Area",
        "LogP",
    ]

    generate_permutations(TEMPLATE, descriptors, pubchem_molecules)
