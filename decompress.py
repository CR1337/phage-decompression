from biotite.sequence.io import fasta, gff
import biotite.sequence.graphics as graphics
from biotite.sequence import NucleotideSequence, Feature, Location, Annotation, AnnotatedSequence, find_subsequence
from typing import Generator, List, Tuple
from dataclasses import dataclass
from itertools import zip_longest
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


IN_FASTA_FILENAME: str = "phiX.fasta"
IN_GFF_FILENAME: str = "phiX.gff"

OUT_FASTA_FILENAME: str = "phiX_decompressed.fasta"
OUT_GFF_FILENAME: str = "phiX_decompressed.gff"

PLOT_FILENAME: str = "plot.png"


@dataclass
class FastaParameters:
    header: str
    chars_per_line: int


@dataclass
class GFFParameters:
    sequence_id: str
    source: str


def load_genome(sequence_filename: str, annotation_filename: str) -> Tuple[AnnotatedSequence, FastaParameters, GFFParameters]:
    fasta_file = fasta.FastaFile.read(sequence_filename)
    genome = fasta.get_sequence(fasta_file, seq_type=NucleotideSequence)

    fasta_parameters = FastaParameters(
        list(fasta_file._entries.keys())[0],
        len(fasta_file.lines[1])
    )

    gff_file = gff.GFFFile.read(annotation_filename)
    annotation = gff.get_annotation(gff_file)

    gff_feature_lines = [line for line in gff_file.lines if not line.startswith("##")]
    sequence_id, source, *_ = gff_feature_lines[0].split("\t")
    gff_parameters = GFFParameters(sequence_id, source)

    return AnnotatedSequence(annotation, genome), fasta_parameters, gff_parameters


def sub_sequences(genome: AnnotatedSequence) -> Generator[Tuple[int, int, Feature], None, None]:
    features = sorted(genome.annotation._features, key=lambda f: next(iter(f.locs)).first)
    
    prev_end = 0
    feature_0_first = next(iter(features[0].locs)).first

    if feature_0_first > 1:
        prev_end = feature_0_first
        gap_location = Location(1, feature_0_first)
        yield 1, feature_0_first, Feature('gap', [gap_location])

    for feature, next_feature in zip_longest(features, features[1:]):
        location = next(iter(feature.locs))

        start, end = location.first, location.last + 1
        prev_end = end
        yield start, end, feature

        if next_feature is not None:
            next_location = next(iter(next_feature.locs))
            next_start = next_location.first
            if next_start > end:
                prev_end = next_start
                gap_location = Location(end, next_start, location.strand)
                yield end, next_start, Feature('gap', [gap_location])

    if prev_end < len(genome.sequence) + 1:
        gap_location = Location(prev_end, len(genome.sequence) + 1)
        yield prev_end, len(genome.sequence) + 1, Feature('gap', [gap_location])


def de_overlap(genome: AnnotatedSequence) -> AnnotatedSequence:
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
    sequences: str | List [str], 
    key: str, 
    *, 
    start: int = 0,
    end: int = -1, 
    first_k: int | None = None, 
    frame: int | None = None
) -> Annotation:
    if isinstance(sequences, str):
        sequences = [sequences]

    features = []

    for s in sequences:
        starts = find_subsequence(genome.sequence[start:end], NucleotideSequence(s))
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


def detect_start_codons(genome: AnnotatedSequence) -> Annotation:
    return detect_feature(genome, ["ATG", "GTG"], "start_codon")


def detect_stop_codons(genome: AnnotatedSequence) -> Annotation:
    return detect_feature(genome, ["TAA", "TAG", "TGA"], "stop_codon")


def detect_rbs(genome: AnnotatedSequence) -> Annotation:
    return detect_feature(genome, ["AGGAGG", "GGAG", "AAGGAG"], "RBS")


def store_genome(genome: AnnotatedSequence, fasta_parameters: FastaParameters, gff_parameters: GFFParameters, sequence_filename: str, annotation_filename: str):
    fasta_file = fasta.FastaFile(fasta_parameters.chars_per_line)
    fasta.set_sequence(fasta_file, genome.sequence, fasta_parameters.header)
    fasta_file.write(sequence_filename)

    gff_file = gff.GFFFile()
    gff.set_annotation(gff_file, genome.annotation, gff_parameters.sequence_id, gff_parameters.source)
    gff_file.write(annotation_filename)


def plot(genome: AnnotatedSequence, decompressed_genome: AnnotatedSequence, filename: str):
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


def main():
    genome, fasta_parameters, gff_parameters = load_genome(IN_FASTA_FILENAME, IN_GFF_FILENAME)
    decompressed_genome = de_overlap(genome)
    
    functional_annotation = detect_start_codons(genome) + detect_stop_codons(genome) + detect_rbs(genome)
    decompressed_functional_features = detect_start_codons(decompressed_genome) + detect_stop_codons(decompressed_genome) + detect_rbs(decompressed_genome)
    
    genome = AnnotatedSequence(genome.annotation + functional_annotation, genome.sequence)
    decompressed_genome = AnnotatedSequence(decompressed_genome.annotation + decompressed_functional_features, decompressed_genome.sequence)


    new_fasta_parameters = FastaParameters(
        f"{fasta_parameters.header} (decompressed)",
        fasta_parameters.chars_per_line
    )
    store_genome(decompressed_genome, new_fasta_parameters, gff_parameters, OUT_FASTA_FILENAME, OUT_GFF_FILENAME)
    plot(genome, decompressed_genome, PLOT_FILENAME)


if __name__ == "__main__":
    main()
