import os
import json
import statistics
import pandas as pd
from tqdm import tqdm
from dataclasses import dataclass
from typing import List
from numbers import Number

import matplotlib.pyplot as plt
import seaborn as sns

from biotite.sequence.io import fasta
from biotite.sequence import NucleotideSequence


PROCESSED_NAMES_FILENAME: str = "processed_names.txt"
DESCRIPTION_FILENAME: str = "analysis.txt"
ANALYSIS_DATA_FILENAME: str = "analysis.csv"
GENOME_DIRECTORY: str = "GenomesDB"


@dataclass
class GenomeAnalysis:
    name: str
    original_length: int
    decompressed_length: int

    @property
    def compression_ratio(self) -> float:
        return self.decompressed_length / self.original_length
    
    @property
    def difference(self) -> int:
        return self.decompressed_length - self.original_length
    
    def save_as_json(self, filename: str):
        data = {
            'name': self.name,
            'original_length': self.original_length,
            'decompressed_length': self.decompressed_length,
            'compression_ratio': self.compression_ratio,
            'difference': self.difference
        }
        with open(filename, 'w') as file:
            json.dump(data, file, indent=4)
    

def analyze(name: str) -> GenomeAnalysis:
    fasta_filename = os.path.join(GENOME_DIRECTORY, f"{name}.fasta")
    decompressed_fasta_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed.fasta")

    fasta_file = fasta.FastaFile.read(fasta_filename)
    decompressed_fasta_file = fasta.FastaFile.read(decompressed_fasta_filename)

    genome = fasta.get_sequence(fasta_file, seq_type=NucleotideSequence)
    decompressed_genome = fasta.get_sequence(decompressed_fasta_file, seq_type=NucleotideSequence)

    result = GenomeAnalysis(name, len(genome), len(decompressed_genome))

    result_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_analysis.json")
    result.save_as_json(result_filename)

    return result


def describe(data: List[Number], title: str) -> str:
    string = ""
    string += title + "\n"
    string += "-" * len(title) + "\n"
    string += f"Mean: {statistics.mean(data)}" + "\n"
    string += f"Median: {statistics.median(data)}" + "\n"
    string += f"StdDev: {statistics.stdev(data)}" + "\n"
    string += f"Min: {min(data)}" + "\n"
    string += f"Max: {max(data)}" + "\n"
    string += "" + "\n"
    return string


def plot(data: List[Number], title: str, ax):
    sns.kdeplot(data, fill=True, color="skyblue", label="KDE", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Value")
    ax.set_ylabel("Density")


def main():
    with open(PROCESSED_NAMES_FILENAME, 'r') as file:
        names = (n.strip() for n in file.readlines())
        names = [n for n in names if len(n) > 0]

    analyses = [
        analyze(name) for name 
        in tqdm(names, total=len(names), desc="Analyzing genomes")
    ]

    names = [a.name for a in analyses]
    original_lengths = [a.original_length for a in analyses]
    decompressed_lengths = [a.decompressed_length for a in analyses]
    compression_ratios = [a.compression_ratio for a in analyses]
    differences = [a.difference for a in analyses]

    description = ""
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))
    axs_flat = axs.ravel()
    for data, title, ax in zip(
        [original_lengths, decompressed_lengths, compression_ratios, differences],
        ["Original lenghts", "Decompressed lengths", "Compression ratios", "Differences"],
        axs_flat
    ):
        description += describe(data, title)
        plot(data, title, ax)

    fig.tight_layout()
    filename = os.path.join("plots", "analysis.png")
    fig.savefig(filename)

    print(description)
    with open(DESCRIPTION_FILENAME, 'w') as file:
        file.write(description)

    df = pd.DataFrame({
        'name': names,
        'original_length': original_lengths,
        'decompressed_length': decompressed_lengths,
        'compression_ratio': compression_ratios,
        'difference': differences
    })
    df.to_csv(ANALYSIS_DATA_FILENAME, index=False)


if __name__ == "__main__":
    main()
