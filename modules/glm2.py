from modules.model import AnnotatedSequence
from modules.decompress import sub_sequences

def create_glm2_string_from_features(sequence: AnnotatedSequence) -> str:
    result = ""

    for _0, _1, feature in sub_sequences(sequence):
        result += "<+>" if feature.location.strand == 1 else "<->"

        if feature.type =='CDS':
            result += str(feature.translate(sequence.sequence, cds=False))

        else:
            result += str(feature.extract(sequence.sequence)).lower()

    return result.replace("*", "")


def create_glm2_string_from_dna(sequence: AnnotatedSequence) -> str:
    return f"<+>{sequence.sequence}".replace("*", "").lower()
