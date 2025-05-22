from dataclasses import dataclass
from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


@dataclass
class AnnotatedSequence:
    record: SeqRecord
    features: List[SeqFeature]

    @property
    def sequence(self) -> Seq:
        return self.record.seq

    def __len__(self) -> int:
        return len(self.record.seq)