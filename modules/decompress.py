from typing import Generator, Tuple
from itertools import zip_longest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from modules.model import AnnotatedSequence



def sub_sequences(genome: AnnotatedSequence) -> Generator[Tuple[int, int, SeqFeature], None, None]:
    """
    Yield non-overlapping subsequences and gaps between features in order.
    """
    features = sorted(genome.features, key=lambda f: f.location.start)

    prev_end = 0
    feature_0_first = features[0].location.start

    # Handle initial gap before first feature
    if feature_0_first > 1:
        prev_end = feature_0_first
        gap_location = SimpleLocation(1, feature_0_first, strand=1)
        yield 1, feature_0_first, SeqFeature(gap_location, type='gap')

    for feature, next_feature in zip_longest(features, features[1:]):
        location = feature.location

        start, end = location.start, location.end + 1
        prev_end = end
        yield start, end, feature

        # Yield a gap if there's space between current and next feature
        if next_feature is not None:
            next_location = next_feature.location
            next_start = next_location.start
            if next_start > end:
                prev_end = next_start
                gap_location = SimpleLocation(end, next_start, strand=1)
                yield end, next_start, SeqFeature(gap_location, type='gap')

    # Handle terminal gap if feature ends before end of sequence
    if prev_end < len(genome.sequence) + 1:
        gap_location = SimpleLocation(prev_end, len(genome.sequence) + 1, strand=1)
        yield prev_end, len(genome.sequence) + 1, SeqFeature(gap_location, 'gap')


def decompress(genome: AnnotatedSequence) -> AnnotatedSequence:
    """
    Decompress overlapping genome features into a new linearized sequence.
    """
    decompressed_sequence = ""
    decompressed_features = []

    for start, end, feature in sub_sequences(genome):
        length = end - start
        decompressed_sequence += str(genome.sequence[start:end])

        location = feature.location
        new_first = len(decompressed_sequence) - length
        new_last = new_first + length - 1
        new_location = SimpleLocation(new_first, new_last, location.strand)
        new_feature = SeqFeature(new_location, type=feature.type, qualifiers=feature.qualifiers)
        decompressed_features.append(new_feature)

    decompressed_genome = SeqRecord(
        Seq(decompressed_sequence), 
        id=genome.record.id,
        name=genome.record.name, 
        description=genome.record.description
    )

    return AnnotatedSequence(decompressed_genome, decompressed_features)


def count_overlapping_cds(sequence: AnnotatedSequence) -> int:
    intervals = sorted(sequence.features, key=lambda f: f.location.start)

    overlapping = set()
    prev_start, prev_end = intervals[0].location.start, intervals[0].location.end
    
    for i in range(1, len(intervals)):
        curr_start, curr_end = intervals[i].location.start, intervals[i].location.end
        
        if curr_start <= prev_end:
            overlapping.add((prev_start, prev_end))
            overlapping.add((curr_start, curr_end))
            prev_end = max(prev_end, curr_end)
        else:
            prev_start, prev_end = curr_start, curr_end
    
    return len(overlapping)
