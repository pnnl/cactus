"""Script to generate benchmark plots."""

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import font_manager


def plot_smiles_length_histogram(df, column_name="smiles", filename: str = "output") -> None:
    """Plot a histogram of the lengths of SMILES strings in a given DataFrame column.

    Parameters
    ----------
        df (pandas.DataFrame): The DataFrame containing SMILES strings.
        column_name (str, optional): The name of the column containing SMILES strings.
    """
    lengths = df[column_name].apply(len)
    num_samples = len(lengths)
    num_bins = int(math.floor(math.sqrt(num_samples)))
    plt.hist(lengths, bins=num_bins, edgecolor="black")
    plt.xlabel("SMILES Length")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of SMILES String Lengths in '{column_name}' Column")
    plt.save_fig(f"{filename}.png", dpi=300)


def create_qualitative_plot(df: pd.DataFrame, modelname: str):
    """Create plot scoring qualitative results."""
    druglike = df[0:100]["answer"].sum()
    pain = df[100:200]["answer"].sum()
    brenk = df[200:300]["answer"].sum()
    bbb = df[300:400]["answer"].sum()
    gi = df[400:500]["answer"].sum()

    plt.style.use("seaborn-v0_8-colorblind")

    font_path = "/Users/mcna892/Library/Fonts/HackNerdFont-Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)

    plt.rcParams["font.family"] = "monospace"
    plt.rcParams["font.monospace"] = prop.get_name()
    # Define categories
    categories = [
        "Druglikeness",
        "PAINS Filter",
        "Brenk Filter",
        "BB Barrier",
        "GI Absorption",
    ]
    incorrect = [100 - val for val in [druglike, pain, brenk, bbb, gi]]
    df_summary = pd.DataFrame(
        list(zip(categories, [druglike, pain, brenk, bbb, gi], incorrect, strict=False)),
        columns=["attribute", "correct", "incorrect"],
    )
    plt.figure(figsize=(10, 6))
    plt.bar(df_summary["attribute"], df_summary["correct"], label="Correct")
    plt.bar(
        df_summary["attribute"],
        df_summary["incorrect"],
        bottom=df_summary["correct"],
        label="Incorrect",
    )
    plt.xlabel("Question Type")
    plt.ylabel("# of Questions")
    plt.xticks(rotation=45, ha="right")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{modelname}_quantitative.png", dpi=300)


def create_quantitative_plot(df: pd.DataFrame, modelname: str):
    """Create plot scoring quantitative results."""
    molwt = df[0:100]["answer"].sum()
    qed = df[100:200]["answer"].sum()
    sa = df[200:300]["answer"].sum()
    tpsa = df[300:400]["answer"].sum()
    logp = df[400:500]["answer"].sum()

    plt.style.use("seaborn-v0_8-colorblind")

    font_path = "/Users/mcna892/Library/Fonts/HackNerdFont-Regular.ttf"  # Your font path goes here
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)

    plt.rcParams["font.family"] = "monospace"
    plt.rcParams["font.monospace"] = prop.get_name()
    categories = ["MolWt", "QED", "SA", "TPSA", "LogP"]
    incorrect = [100 - val for val in [molwt, qed, sa, tpsa, logp]]
    df_summary = pd.DataFrame(
        list(zip(categories, [molwt, qed, sa, tpsa, logp], incorrect, strict=False)),
        columns=["attribute", "correct", "incorrect"],
    )

    plt.figure(figsize=(10, 6))

    plt.bar(df_summary["attribute"], df_summary["correct"], label="Correct")
    plt.bar(
        df_summary["attribute"],
        df_summary["incorrect"],
        bottom=df_summary["correct"],
        label="Incorrect",
    )
    plt.xlabel("Question Type")
    plt.ylabel("# of Questions")
    plt.xticks(rotation=45, ha="right")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{modelname}_quantitative.png", dpi=300)


def create_combined_plot(df: pd.DataFrame, model_name: str):
    """Create plot scoring quantitative and qualitative results."""
    druglike = df[0:100]["answer"].sum()
    pain = df[100:200]["answer"].sum()
    brenk = df[200:300]["answer"].sum()
    bbb = df[300:400]["answer"].sum()
    gi = df[400:500]["answer"].sum()
    molwt = df[500:600]["answer"].sum()
    qed = df[600:700]["answer"].sum()
    sa = df[700:800]["answer"].sum()
    tpsa = df[800:900]["answer"].sum()
    logp = df[900:1000]["answer"].sum()

    plt.style.use("seaborn-v0_8-colorblind")

    font_path = "/Users/mcna892/Library/Fonts/HackNerdFont-Regular.ttf"  # Your font path goes here
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)

    plt.rcParams["font.family"] = "monospace"
    plt.rcParams["font.monospace"] = prop.get_name()
    # Define categories
    categories = [
        "Druglikeness",
        "PAINS Filter",
        "Brenk Filter",
        "BB Barrier",
        "GI Absorption",
        "MolWt",
        "QED",
        "SA",
        "TPSA",
        "LogP",
    ]
    incorrect = [100 - val for val in [druglike, pain, brenk, bbb, gi, molwt, qed, sa, tpsa, logp]]
    df_summary = pd.DataFrame(
        list(
            zip(
                categories,
                [druglike, pain, brenk, bbb, gi, molwt, qed, sa, tpsa, logp],
                incorrect,
                strict=False,
            )
        ),
        columns=["attribute", "correct", "incorrect"],
    )
    mid_index = (len(categories) - 1) // 2
    x = np.arange(len(categories))
    bar_width = 1
    plt.figure(figsize=(10, 6))
    plt.axvline(x[mid_index] + bar_width / 2, color="black", linestyle="--", linewidth=1)
    plt.text(2, 107, "Qualitative", ha="center", va="center", fontsize=15)
    plt.text(7, 107, "Quantitative", ha="center", va="center", fontsize=15)

    plt.bar(df_summary["attribute"], df_summary["correct"], label="Correct")
    plt.bar(
        df_summary["attribute"],
        df_summary["incorrect"],
        bottom=df_summary["correct"],
        label="Incorrect",
    )
    plt.xlabel("Question Type")
    plt.ylabel("# of Questions")
    plt.xticks(rotation=45, ha="right")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{model_name}_Combined.png", dpi=300)


def create_aggregate_plot(dfs: list, model_names: list):
    """Create plot for comparing correct results from each model per category."""
    categories = [
        "Druglikeness",
        "PAINS Filter",
        "Brenk Filter",
        "BB Barrier",
        "GI Absorption",
        "MolWt",
        "QED",
        "SA",
        "TPSA",
        "LogP",
    ]

    plt.style.use("seaborn-v0_8-colorblind")
    font_path = "/Users/mcna892/Library/Fonts/HackNerdFont-Regular.ttf"  # Your font path goes here
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams["font.family"] = "monospace"
    plt.rcParams["font.monospace"] = prop.get_name()

    plt.figure(figsize=(10, 6))

    num_models = len(model_names)
    bar_width = 0.7 / num_models
    x = np.arange(len(categories))

    colors = sns.color_palette("colorblind", 10)

    for i, (df, model_name) in enumerate(zip(dfs, model_names, strict=False)):
        correct = [df[j * 100 : (j + 1) * 100]["answer"].sum() for j in range(10)]
        plt.bar(
            x - 0.4 + i * bar_width,
            correct,
            bar_width,
            color=colors[i],
            edgecolor="black",
            label=f"{model_name}",
        )

    plt.xlabel("Question Type")
    plt.ylabel("# of Questions")
    plt.xticks(x, categories, rotation=45, ha="right")
    plt.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.05),
        ncol=4,
    )
    plt.tight_layout()
    plt.savefig("Combined.png", dpi=600)
