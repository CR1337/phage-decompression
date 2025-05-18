import json
import statistics
from numbers import Number
from dataclasses import dataclass
from typing import List

import pandas as pd
from biotite.sequence import NucleotideSequence


@dataclass
class GenomeAnalysis:
    name: str
    original_length: int
    decompressed_length: int

    @property
    def compression_ratio(self) -> float:
        return self.decompressed_length / self.original_length
    
    @property
    def compression_factor(self) -> float:
        return self.original_length / self.decompressed_length
    
    @property
    def difference(self) -> int:
        return self.decompressed_length - self.original_length
    
    def save_as_json(self, filename: str):
        data = {
            'name': self.name,
            'original_length': self.original_length,
            'decompressed_length': self.decompressed_length,
            'compression_ratio': self.compression_ratio,
            'compression_factor': self.compression_factor,
            'difference': self.difference
        }
        with open(filename, 'w') as file:
            json.dump(data, file, indent=4)
    

def analyze(name: str, genome: NucleotideSequence, decompressed_genome: NucleotideSequence) -> GenomeAnalysis:
    return GenomeAnalysis(name, len(genome), len(decompressed_genome))


def _value_summary(data: List[Number], title: str) -> str:
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


def summary(analyses: List[GenomeAnalysis]) -> str:
    values = [
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
    result = ""
    for data, title in zip(values, titles):
        result += _value_summary(data, title)
    return result


def make_dataframe(analyses: List[GenomeAnalysis]) -> pd.DataFrame:
    names = [a.name for a in analyses]
    original_lengths = [a.original_length for a in analyses]
    decompressed_lengths = [a.decompressed_length for a in analyses]
    compression_factors = [a.compression_factor for a in analyses]
    differences = [a.difference for a in analyses]

    df = pd.DataFrame({
        'name': names,
        'original_length': original_lengths,
        'decompressed_length': decompressed_lengths,
        'compression_factor': compression_factors,
        'difference': differences
    })

    return df
