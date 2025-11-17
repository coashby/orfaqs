"""
Proteins
"""

import logging

from orfaqs.lib.core.aminoacids import (
    AminoAcid,
    AminoAcidSequence,
)

_logger = logging.getLogger(__name__)


class Protein(AminoAcidSequence):
    """Protein"""

    def __init__(
        self,
        sequence: (AminoAcid | list[AminoAcid] | str) = None,
        name: str = None,
        sequence_complete: bool = False,
    ):
        super().__init__(sequence=sequence, name=name)

        self._sequence_complete = sequence_complete
        if isinstance(sequence, Protein):
            self._sequence_complete = sequence.sequence_complete

    @property
    def sequence_complete(self) -> bool:
        return self._sequence_complete

    def add_amino_acid(self, amino_acid: (AminoAcid | str)):
        if self._sequence_complete:
            return
        if isinstance(amino_acid, AminoAcid):
            amino_acid = amino_acid.symbol

        self._sequence_str += amino_acid

    def end_sequence(self):
        self._sequence_complete = True
