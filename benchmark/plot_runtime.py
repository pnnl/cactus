"""Script to create the Accuracy vs Time Plot."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

color_palette = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#8C510A",
    "#56B4E9",
    "#E69F00",
    "#6A3D9A",
    "#808080",
    "#000000",
]


def plot_accuracy_v_time_old(file_path):
    """Load in data to plot accuracy vs time."""
    df = pd.read_csv(file_path)

    categories = ["Prompt", "Node", "Question Set"]

    _, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, category in enumerate(categories):
        sns.scatterplot(
            x="Time",
            y="Accuracy",
            hue="Model",
            data=df,
            style=category,
            palette=color_palette,
            ax=axes[i],
            s=200,
        )

        axes[i].set_title(f"Style by{category}")

    plt.tight_layout()

    plt.show()


def plot_accuracy_v_time_new(file_path):
    """Plot the accuracy vs time under different conditions."""
    df = pd.read_csv(file_path)

    g = sns.FacetGrid(
        df,
        col="Question Set",
        col_order=["Qualitative", "Quantitative", "Full"],
        row="Prompt",
        row_order=["Minimal", "Domain"],
        sharex=False,
        height=4,
        aspect=1.2,
    )

    g.map(
        sns.scatterplot,
        "Time",
        "Accuracy",
        data=df,
        hue="Model",
        style="Device",
        palette=color_palette,
        s=100,
    )
    g.add_legend()
    g.tight_layout()
    plt.savefig("./benchmark_files/acc_vs_time_v3.png", dpi=300)


if __name__ == "__main__":
    plot_accuracy_v_time_new("./benchmark_files/runtime.csv")
