from math import sqrt, ceil
from numbers import Number
from typing import Dict, List, Tuple, Callable

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from biotite.sequence import AnnotatedSequence, Feature
import biotite.sequence.graphics as graphics

from modules.analysis import GenomeAnalysis


def plot_decompressed_plasmid(
    genome: AnnotatedSequence, 
    decompressed_genome: AnnotatedSequence, 
    filename: str,
    feature_format: Dict[str, Tuple[bool, str, str, Callable[[Feature], str]]]
):
    def feature_formatter(feature: Feature) -> Tuple[bool, str, str, str | None]:
        f = feature_format[feature.key]
        return f[0], f[1], f[2], f[3](feature)

    fig, axs = plt.subplots(
        1, 2, subplot_kw={'projection': 'polar'}, figsize=(16, 8)
    )
    
    graphics.plot_plasmid_map(
        axs[0], 
        genome.annotation, 
        len(genome.sequence), 
        label="compressed",
        feature_formatter=feature_formatter
    )

    graphics.plot_plasmid_map(
        axs[1], 
        decompressed_genome.annotation, 
        len(decompressed_genome.sequence), 
        label="decompressed",
        feature_formatter=feature_formatter
    )

    legend_elements = [
        Patch(facecolor=v[1], edgecolor='black', label=k) 
        for k, v in feature_format.items()
    ]
    fig.legend(
        handles=legend_elements, loc='lower center', title='Feature Types'
    )
    fig.savefig(filename)


def plot_analysis(analyses: List[GenomeAnalysis], filename: str, log_scale: bool=False):
    data_list = [
        [a.original_length for a in analyses],
        [a.decompressed_length for a in analyses],
        [a.compression_factor for a in analyses],
        [a.difference for a in analyses]
    ]
    titles = [
        "Original lenghts", 
        "Decompressed lengths", 
        "Compression factors", 
        "Differences"
    ]
    if log_scale:
        data_list = [np.log10(data) for data in data_list]
    kde_plots(data_list, titles, filename)


def kde_plots(data_list: List[List[Number]], titles: List[str], filename: str):
    width = height = ceil(sqrt(len(data_list)))
    fig, axs = plt.subplots(width, height, figsize=(width * 6, height * 6))
    axs_flat = axs.ravel()
    for data, title, ax in zip(data_list, titles, axs_flat):
        sns.kdeplot(data, fill=True, color='skyblue', ax=ax)
        ax.set_title(title)
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")
    fig.tight_layout()
    fig.savefig(filename)
