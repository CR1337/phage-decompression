import os
from dataclasses import dataclass
from typing import Tuple

from biotite.sequence import AnnotatedSequence, NucleotideSequence
from biotite.sequence.io import fasta, gff, genbank


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


@dataclass
class SequenceFile:
    name: str
    directory: str

    @property
    def fasta_filename(self) -> str:
        return os.path.join(self.directory, f"{self.name}.fasta")
    
    @property
    def decompressed_fasta_filename(self) -> str:
        return os.path.join(self.directory, self.name, f"{self.name}_decompressed.fasta")
    
    @property 
    def gff_filename(self) -> str:
        return os.path.join(self.directory, self.name, f"{self.name}.gff")
    
    @property
    def decompressed_gff_filename(self) -> str:
        return os.path.join(self.directory, self.name, f"{self.name}_decompressed.gff")
    
    @property
    def has_fasta_file(self) -> bool:
        return os.path.exists(self.fasta_filename)
    
    @property
    def has_gff_file(self) -> bool:
        return os.path.exists(self.gff_filename)
    
    @property
    def has_decompressed_fasta_file(self) -> bool:
        return os.path.exists(self.decompressed_fasta_filename)
    
    @property
    def has_decompressed_gff_file(self) -> bool:
        return os.path.exists(self.decompressed_gff_filename)
    
    def make_path(self, filename: str) -> str:
        return os.path.join(self.directory, self.name, filename)


def load_annotated_sequence(genbank_filename: str) -> Tuple[AnnotatedSequence, FastaParameters, GFFParameters]:
    """
    Load a genome sequence and its annotation from FASTA and GFF files.
    """
    genbank_file = genbank.GenBankFile.read(genbank_filename)
    return genbank.get_annotated_sequence(genbank_file), FastaParameters(genbank.get_definition(genbank_file), 80), GFFParameters("NO_LOCUS", genbank.get_source(genbank_file))


    fasta_file = fasta.FastaFile.read(fasta_filename)
    genome = fasta.get_sequence(fasta_file, seq_type=NucleotideSequence)

    # Extract FASTA parameters
    fasta_parameters = FastaParameters(
        list(fasta_file._entries.keys())[0],
        len(fasta_file.lines[1])
    )

    gff_file = gff.GFFFile.read(gff_filename)
    annotation = gff.get_annotation(gff_file)

    # Parse GFF metadata from first non-comment line
    gff_feature_lines = [line for line in gff_file.lines if not line.startswith("##")]
    sequence_id, source, *_ = gff_feature_lines[0].split("\t")
    gff_parameters = GFFParameters(sequence_id, source)

    return AnnotatedSequence(annotation, genome), fasta_parameters, gff_parameters


def load_sequence(file: SequenceFile, decompressed: bool = False) -> NucleotideSequence:
    fasta_file = fasta.FastaFile.read(file.decompressed_fasta_filename if decompressed else file.fasta_filename)
    genome = fasta.get_sequence(fasta_file, seq_type=NucleotideSequence)
    return genome


def store_decompressed_sequence(genome: AnnotatedSequence, fasta_filename: str, gff_filename: str, fasta_parameters: FastaParameters, gff_parameters: GFFParameters):
    fasta_file = fasta.FastaFile(fasta_parameters.chars_per_line)
    fasta.set_sequence(fasta_file, genome.sequence, fasta_parameters.header)
    fasta_file.write(fasta_filename)

    gff_file = gff.GFFFile()
    gff.set_annotation(gff_file, genome.annotation, gff_parameters.sequence_id, gff_parameters.source)
    gff_file.write(gff_filename)
