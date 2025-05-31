from dataclasses import dataclass
from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


@dataclass
class AnnotatedSequence:
    """
    A container for a biological sequence along with its associated features.

    Attributes:
        record (SeqRecord): The biological sequence record.
        features (List[SeqFeature]): A list of annotated features (e.g., CDS, gene).
    """

    record: SeqRecord
    features: List[SeqFeature]

    @property
    def sequence(self) -> Seq:
        """
        Get the nucleotide or amino acid sequence of the record.

        Returns:
            Seq: The sequence from the underlying SeqRecord.
        """
        return self.record.seq

    def __len__(self) -> int:
        """
        Get the length of the sequence.

        Returns:
            int: The number of bases or residues in the sequence.
        """
        return len(self.record.seq)
