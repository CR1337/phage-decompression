import os
from math import sqrt, ceil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Seq import Seq


def compute_similarity(sequence_1: Seq, sequence_2: Seq, k: int = 10) -> float:
    """
    Compute the Jaccard similarity between two sequences using k-mer sets.

    The similarity is defined as the size of the intersection divided by
    the size of the union of k-mers from both sequences.

    Parameters:
        sequence_1 (Seq): The first sequence.
        sequence_2 (Seq): The second sequence.

    Returns:
        float: Jaccard similarity score between the two sequences.
    """
    kmers1 = {sequence_1[i : i + k] for i in range(len(sequence_1) - k + 1)}
    kmers2 = {sequence_2[i : i + k] for i in range(len(sequence_2) - k + 1)}
    intersection = kmers1 & kmers2
    union = kmers1 | kmers2
    similarity = len(intersection) / len(union) if union else 0.0
    return similarity


def plot_similarities(df: pd.DataFrame, filename: str, cluster: bool):
    """
    Plot a similarity matrix from a DataFrame and save it as an image.

    The DataFrame must contain the columns: 'filename_1', 'filename_2', and 'similarity'.
    Sequence names are derived from file basenames. Optionally applies hierarchical clustering.

    Parameters:
        df (pd.DataFrame): DataFrame containing pairwise similarity data.
        filename (str): Path to save the output image (e.g., heatmap or cluster map).
        cluster (bool): If True, uses seaborn clustermap; otherwise, uses a regular heatmap.

    Returns:
        None
    """
    df["sequence 1"] = df["filename_1"].apply(
        lambda x: os.path.splitext(os.path.basename(x))[0]
    )
    df["sequence 2"] = df["filename_2"].apply(
        lambda x: os.path.splitext(os.path.basename(x))[0]
    )
    df.loc[df["filename_1"] == df["filename_2"], "similarity"] = 0.0
    sim_matrix = df.pivot(index="sequence 1", columns="sequence 2", values="similarity")
    sim_matrix = sim_matrix.fillna(0)
    size = min(max(ceil(sqrt((len(df) // 3))), 10) + 2, 16_000)
    plt.figure(figsize=(size, size))
    if cluster:
        sns.clustermap(sim_matrix, cmap="viridis", figsize=(size, size))
    else:
        sns.heatmap(
            sim_matrix,
            cmap="viridis",
            annot=False,
            square=True,
            xticklabels=False,
            yticklabels=False,
        )
    plt.title("Sequence Similarity Matrix")
    plt.tight_layout()
    plt.savefig(filename)
