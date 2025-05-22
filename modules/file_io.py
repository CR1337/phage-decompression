from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

from modules.model import AnnotatedSequence


def load_annotated_sequence_from_fasta(fasta_filename: str, gff_filename: str) -> AnnotatedSequence:
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))
    sequence = list(sequences.values())[0]
    with open(gff_filename, 'r') as gff_file:
        features = list(next(GFF.parse(gff_file, base_dict=sequences)).features)
    return AnnotatedSequence(sequence, features)


def load_annotated_sequence_from_genbank(genbank_filename: str) -> AnnotatedSequence:
    record = next(SeqIO.parse(genbank_filename, 'genbank'))
    sequence = record
    features = [f for f in record.features if f.type != 'source']
    return AnnotatedSequence(sequence, features)


def load_sequence_from_fasta(fasta_filename: str) -> SeqRecord:
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))
    sequence = list(sequences.values())[0]
    return sequence


def load_sequence_from_genbank(genbank_filename: str) -> SeqRecord:
    record = next(SeqIO.parse(genbank_filename, 'genbank'))
    return record


def store_annotated_sequence_as_fasta(sequence: AnnotatedSequence, fasta_filename: str, gff_filename: str):
    record = SeqRecord(
        seq=sequence.sequence,
        id=sequence.record.id,
        name=sequence.record.name,
        description=sequence.record.description,
        annotations=sequence.record.annotations,
        features=sequence.features
    )

    with open(fasta_filename, "w") as fasta_file:
        SeqIO.write(record, fasta_file, "fasta")

    with open(gff_filename, "w") as gff_file:
        GFF.write([record], gff_file)
