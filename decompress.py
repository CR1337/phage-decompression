import os
import json
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Generator, List, Tuple
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from tqdm import tqdm

from biotite.sequence import (
    NucleotideSequence, Feature, Location, Annotation,
    AnnotatedSequence, find_subsequence, ProteinSequence
)
from biotite.sequence.io import fasta, gff
import biotite.sequence.graphics as graphics

import warnings

# Ignore specific UserWarning about FASTA data in GFF files from Biotite
warnings.filterwarnings("ignore", message="Biotite does not support FASTA data mixed into GFF files")



PROCESSED_NAMES_FILENAME: str = "processed_names.txt"
GENOME_DIRECTORY: str = "GenomesDB"


# Holds parameters for FASTA output
@dataclass
class FastaParameters:
    header: str
    chars_per_line: int


# Holds parameters for GFF output
@dataclass
class GFFParameters:
    sequence_id: str
    source: str


def load_genome(sequence_filename: str, annotation_filename: str) -> Tuple[AnnotatedSequence, FastaParameters, GFFParameters]:
    """
    Load a genome sequence and its annotation from FASTA and GFF files.
    """
    fasta_file = fasta.FastaFile.read(sequence_filename)
    genome = fasta.get_sequence(fasta_file, seq_type=NucleotideSequence)

    # Extract FASTA parameters
    fasta_parameters = FastaParameters(
        list(fasta_file._entries.keys())[0],
        len(fasta_file.lines[1])
    )

    gff_file = gff.GFFFile.read(annotation_filename)
    annotation = gff.get_annotation(gff_file)

    # Parse GFF metadata from first non-comment line
    gff_feature_lines = [line for line in gff_file.lines if not line.startswith("##")]
    sequence_id, source, *_ = gff_feature_lines[0].split("\t")
    gff_parameters = GFFParameters(sequence_id, source)

    return AnnotatedSequence(annotation, genome), fasta_parameters, gff_parameters


def sub_sequences(genome: AnnotatedSequence) -> Generator[Tuple[int, int, Feature], None, None]:
    """
    Yield non-overlapping subsequences and gaps between features in order.
    """
    features = sorted(genome.annotation._features, key=lambda f: next(iter(f.locs)).first)
    
    prev_end = 0
    feature_0_first = next(iter(features[0].locs)).first

    # Handle initial gap before first feature
    if feature_0_first > 1:
        prev_end = feature_0_first
        gap_location = Location(1, feature_0_first)
        yield 1, feature_0_first, Feature('gap', [gap_location])

    for feature, next_feature in zip_longest(features, features[1:]):
        location = next(iter(feature.locs))

        start, end = location.first, location.last + 1
        prev_end = end
        yield start, end, feature

        # Yield a gap if there's space between current and next feature
        if next_feature is not None:
            next_location = next(iter(next_feature.locs))
            next_start = next_location.first
            if next_start > end:
                prev_end = next_start
                gap_location = Location(end, next_start, location.strand)
                yield end, next_start, Feature('gap', [gap_location])

    # Handle terminal gap if feature ends before end of sequence
    if prev_end < len(genome.sequence) + 1:
        gap_location = Location(prev_end, len(genome.sequence) + 1)
        yield prev_end, len(genome.sequence) + 1, Feature('gap', [gap_location])


def de_overlap(genome: AnnotatedSequence) -> AnnotatedSequence:
    """
    Decompress overlapping genome features into a new linearized sequence.
    """
    decompressed_sequence = ""
    decompressed_features = []

    for start, end, feature in sub_sequences(genome):
        length = end - start
        decompressed_sequence += str(genome.sequence[start:end])

        location = next(iter(feature.locs))
        new_first = len(decompressed_sequence) - length
        new_last = new_first + length - 1
        new_location = Location(new_first, new_last, location.strand, location.defect)
        new_feature = Feature(feature.key, [new_location], feature.qual)
        decompressed_features.append(new_feature)

    decompressed_genome = NucleotideSequence(decompressed_sequence)
    decompressed_annotation = Annotation(decompressed_features)

    return AnnotatedSequence(decompressed_annotation, decompressed_genome)


def detect_feature(
    genome: AnnotatedSequence, 
    sequences: str | List[str], 
    key: str, 
    *, 
    start: int = 0,
    end: int = -1, 
    first_k: int | None = None, 
    frame: int | None = None
) -> Annotation:
    """
    Detect specific motifs (e.g., codons, RBS) in the genome and return as features.
    """
    if isinstance(sequences, str):
        sequences = [sequences]

    features = []

    for s in sequences:
        starts = find_subsequence(genome.sequence[start:end], NucleotideSequence(s))

        # Apply frame constraint if specified
        if frame is not None:
            frame %= 3
            starts = [start for start in starts if start % frame == 0]

        for start in starts:
            location = Location(start, start + len(s) - 1)
            feature = Feature(key, [location])
            features.append(feature)

    features = sorted(features, key=lambda f: next(iter(f.locs)).first)
    
    if first_k:
        features = features[:first_k]

    return Annotation(features)


# Feature detection wrappers for specific motifs
def detect_start_codons(genome: AnnotatedSequence) -> Annotation:
    return detect_feature(genome, ["ATG", "GTG"], "start_codon")


def detect_stop_codons(genome: AnnotatedSequence) -> Annotation:
    return detect_feature(genome, ["TAA", "TAG", "TGA"], "stop_codon")


def detect_rbs(genome: AnnotatedSequence) -> Annotation:
    return detect_feature(genome, ["AGGAGG", "GGAG", "AAGGAG"], "RBS")


def store_genome(genome: AnnotatedSequence, fasta_parameters: FastaParameters, gff_parameters: GFFParameters, sequence_filename: str, annotation_filename: str):
    """
    Write the annotated genome sequence to FASTA and GFF files.
    """
    fasta_file = fasta.FastaFile(fasta_parameters.chars_per_line)
    fasta.set_sequence(fasta_file, genome.sequence, fasta_parameters.header)
    fasta_file.write(sequence_filename)

    gff_file = gff.GFFFile()
    gff.set_annotation(gff_file, genome.annotation, gff_parameters.sequence_id, gff_parameters.source)
    gff_file.write(annotation_filename)


def plot_plasmids(genome: AnnotatedSequence, decompressed_genome: AnnotatedSequence, filename: str):
    """
    Plot the plasmid maps of compressed and decompressed genomes.
    """
    def feature_formatter(feature: Feature):
        match feature.key:
            case "CDS":
                return True, "orange", 'black', feature.qual.get('function')
            case "RBS":
                return False, "blue", 'black', None
            case "start_codon":
                return False, "green", 'black', None
            case "stop_codon":
                return False, "red", 'black', None
            case "gap":
                return False, "gray", 'black', None
            case _:
                raise ValueError("Invalid Feature")

    fig, axs = plt.subplots(1, 2, subplot_kw={'projection': 'polar'}, figsize=(16, 8))
    
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
        Patch(facecolor='orange', edgecolor='black', label='CDS'),
        Patch(facecolor='blue', edgecolor='black', label='RBS'),
        Patch(facecolor='green', edgecolor='black', label='Start Codon'),
        Patch(facecolor='red', edgecolor='black', label='Stop Codon'),
        Patch(facecolor='gray', edgecolor='black', label='Gap'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', title='Feature Types')
    fig.savefig(filename)


def translate(genome: AnnotatedSequence) -> List[List[Tuple[ProteinSequence, int, int]]]:
    """
    TODO
    """
    cds_features = (f for f in genome.annotation._features if f.key == "CDS")
    protein_sequences = []
    for feature in cds_features:
        location = next(iter(feature.locs))
        sub_sequence = genome.sequence[location.first:location.last+1]
        aa_sequences, positions = sub_sequence.translate(complete=False)
        protein_sequences.append(sorted(
            [
                (sequence, int(first), int(last)) 
                for sequence, (first, last) in zip(aa_sequences, positions)
            ],
            key=lambda x: len(x[0]),
            reverse=True
        ))
    return protein_sequences


def store_protein_sequences(protein_sequences: List[List[Tuple[ProteinSequence, int, int]]], filename: str):
    """
    TODO
    """
    sequences = [
        [
            {'first': first, 'last': last, 'sequence': str(sequence)[:-1]} 
            for sequence, first, last in sequence_list
        ]
        for sequence_list in protein_sequences
    ]
    with open(filename, 'w') as file:
        json.dump(sequences, file, indent=4)


def process_genome(
    name: str, 
    plot_plasmids: bool = False, translate: bool = False, do_functional_annotation: bool = False) -> bool:
    """
    Full pipeline:
    1. Load compressed genome
    2. Decompress overlapping features
    3. Detect motifs (start/stop codons, RBS)
    4. Save decompressed genome
    5. Generate visualization
    6. Translate into protein sequences
    7. Save protein sequences
    """

    in_fasta_filename = os.path.join(GENOME_DIRECTORY, f"{name}.fasta")
    in_gff_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}.gff")
    out_fasta_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed.fasta")
    out_gff_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed.gff")

    plasmid_plot_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed_plasmids.png")
    out_protein_filename = os.path.join(GENOME_DIRECTORY, name, f"{name}_decompressed_proteins.json")

    if not os.path.exists(in_gff_filename):
        return False

    genome, fasta_parameters, gff_parameters = load_genome(in_fasta_filename, in_gff_filename)
    decompressed_genome = de_overlap(genome)
    
    if do_functional_annotation:
        # Detect additional functional features
        functional_annotation = detect_start_codons(genome) + detect_stop_codons(genome) + detect_rbs(genome)
        decompressed_functional_features = detect_start_codons(decompressed_genome) + detect_stop_codons(decompressed_genome) + detect_rbs(decompressed_genome)
        annotations = genome.annotation + functional_annotation
        decompressed_annotations = decompressed_genome.annotation + decompressed_functional_features
    else:
        annotations = genome.annotation
        decompressed_annotations = decompressed_genome.annotation
    
    # Combine annotations
    genome = AnnotatedSequence(annotations, genome.sequence)
    decompressed_genome = AnnotatedSequence(decompressed_annotations, decompressed_genome.sequence)

    # Save decompressed genome and annotation
    new_fasta_parameters = FastaParameters(
        f"{fasta_parameters.header} (decompressed)",
        fasta_parameters.chars_per_line
    )
    store_genome(decompressed_genome, new_fasta_parameters, gff_parameters, out_fasta_filename, out_gff_filename)

    if plot_plasmids:
        # Generate comparative plot
        plot_plasmids(genome, decompressed_genome, plasmid_plot_filename)

    if translate:
        # Translate decompressed genome
        protein_sequences = translate(decompressed_genome)

        # Save protein sequences
        store_protein_sequences(protein_sequences, out_protein_filename)

    return True


def main():
    directory = Path(GENOME_DIRECTORY)
    names = [
        ".".join(f.name.split(".")[:-1]) 
        for f in directory.rglob('*.fasta')
    ]

    with Pool(processes=cpu_count() - 1) as pool:
        results = list(tqdm(pool.imap(process_genome, names), total=len(names), desc="Processing genomes"))

    with open(PROCESSED_NAMES_FILENAME, 'w') as file:
        for name, result in zip(names, results):
            if result:
                file.write(f"{name}\n")

    print("DONE")


if __name__ == "__main__":
    main()
