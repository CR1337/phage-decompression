import os
from math import sqrt, ceil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Seq import Seq


def compute_similarity(sequence_1: Seq, sequence_2: Seq) -> float:
    k = 10
    kmers1 = {sequence_1[i:i+k] for i in range(len(sequence_1) - k + 1)}
    kmers2 = {sequence_2[i:i+k] for i in range(len(sequence_2) - k + 1)}
    intersection = kmers1 & kmers2
    union = kmers1 | kmers2
    similarity = len(intersection) / len(union) if union else 0.0
    return similarity


def plot_similarities(df: pd.DataFrame, filename: str, cluster: bool):
    df["sequence 1"] = df["filename_1"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df["sequence 2"] = df["filename_2"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df.loc[df["filename_1"] == df["filename_2"], "similarity"] = 0.0
    sim_matrix = df.pivot(index="sequence 1", columns="sequence 2", values="similarity")
    sim_matrix = sim_matrix.fillna(0)
    size = min(max(ceil(sqrt((len(df) // 3))), 10) + 2, 16_000)
    plt.figure(figsize=(size, size))
    if cluster:
        sns.clustermap(sim_matrix, cmap="viridis", figsize=(size, size))
    else:
        sns.heatmap(sim_matrix, cmap="viridis", annot=False, square=True, xticklabels=False, yticklabels=False)
    plt.title("Sequence Similarity Matrix")
    plt.tight_layout()
    plt.savefig(filename)