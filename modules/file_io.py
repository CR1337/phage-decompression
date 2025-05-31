from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

from modules.model import AnnotatedSequence


def load_annotated_sequence_from_fasta(
    fasta_filename: str, gff_filename: str
) -> AnnotatedSequence:
    """
    Load an annotated sequence from a FASTA file and corresponding GFF annotation file.

    Parameters:
        fasta_filename (str): Path to the FASTA file.
        gff_filename (str): Path to the GFF file containing feature annotations.

    Returns:
        AnnotatedSequence: A sequence with associated features parsed from the GFF.
    """
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))
    sequence = list(sequences.values())[0]
    with open(gff_filename, "r") as gff_file:
        features = list(next(GFF.parse(gff_file, base_dict=sequences)).features)
    return AnnotatedSequence(sequence, features)


def load_annotated_sequence_from_genbank(genbank_filename: str) -> AnnotatedSequence:
    """
    Load an annotated sequence from a GenBank file.

    Parameters:
        genbank_filename (str): Path to the GenBank file.

    Returns:
        AnnotatedSequence: A sequence with associated features.
    """
    record = next(SeqIO.parse(genbank_filename, "genbank"))
    sequence = record
    features = [f for f in record.features if f.type != "source"]
    return AnnotatedSequence(sequence, features)


def load_sequence_from_fasta(fasta_filename: str) -> SeqRecord:
    """
    Load a plain sequence (without annotations) from a FASTA file.

    Parameters:
        fasta_filename (str): Path to the FASTA file.

    Returns:
        SeqRecord: The loaded sequence.
    """
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))
    sequence = list(sequences.values())[0]
    return sequence


def load_sequence_from_genbank(genbank_filename: str) -> SeqRecord:
    """
    Load a plain sequence from a GenBank file.

    Parameters:
        genbank_filename (str): Path to the GenBank file.

    Returns:
        SeqRecord: The loaded sequence.
    """
    record = next(SeqIO.parse(genbank_filename, "genbank"))
    return record


def store_annotated_sequence_as_fasta(
    sequence: AnnotatedSequence, fasta_filename: str, gff_filename: str
):
    """
    Store an annotated sequence as a FASTA file and corresponding GFF file.

    Parameters:
        sequence (AnnotatedSequence): The sequence and features to store.
        fasta_filename (str): Output path for the FASTA file.
        gff_filename (str): Output path for the GFF annotation file.
    """
    record = SeqRecord(
        seq=sequence.sequence,
        id=sequence.record.id,
        name=sequence.record.name,
        description=sequence.record.description,
        annotations=sequence.record.annotations,
        features=sequence.features,
    )

    with open(fasta_filename, "w") as fasta_file:
        SeqIO.write(record, fasta_file, "fasta")

    with open(gff_filename, "w") as gff_file:
        GFF.write([record], gff_file)
