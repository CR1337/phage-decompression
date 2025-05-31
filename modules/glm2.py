from modules.model import AnnotatedSequence
from modules.decompress import sub_sequences


def create_glm2_string_from_features(sequence: AnnotatedSequence) -> str:
    """
    Generate a glm2-compatible string representation of a sequence from its features.

    This function processes each feature in the sequence using `sub_sequences`,
    appends a directional tag (`<+>` or `<->`), and includes either:
      - the translated amino acid sequence for CDS features (without stop codons), or
      - the lowercased nucleotide sequence for non-CDS features.

    Parameters:
        sequence (AnnotatedSequence): The annotated sequence to process.

    Returns:
        str: A glm2-formatted string representing the feature-based sequence content.
    """
    result = ""

    for _0, _1, feature in sub_sequences(sequence):
        result += "<+>" if feature.location.strand == 1 else "<->"

        if feature.type == "CDS":
            result += str(feature.translate(sequence.sequence, cds=False))

        else:
            result += str(feature.extract(sequence.sequence)).lower()

    return result.replace("*", "")


def create_glm2_string_from_dna(sequence: AnnotatedSequence) -> str:
    """
    Generate a glm2-compatible string representation of a raw DNA sequence.

    The sequence is prefixed with a `<+>` direction marker, lowercased, and
    any stop codons (`*`) are removed.

    Parameters:
        sequence (AnnotatedSequence): The raw sequence to convert.

    Returns:
        str: A glm2-formatted string representing the unannotated DNA content.
    """
    return f"<+>{sequence.sequence}".replace("*", "").lower()
