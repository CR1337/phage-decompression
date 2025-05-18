from itertools import zip_longest
from typing import Generator, Tuple

from biotite.sequence import (
    NucleotideSequence, Feature, Location, Annotation,
    AnnotatedSequence
)


def _sub_sequences(genome: AnnotatedSequence) -> Generator[Tuple[int, int, Feature], None, None]:
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


def decompress(genome: AnnotatedSequence) -> AnnotatedSequence:
    """
    Decompress overlapping genome features into a new linearized sequence.
    """
    decompressed_sequence = ""
    decompressed_features = []

    for start, end, feature in _sub_sequences(genome):
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
