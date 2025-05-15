"""
Proteins
"""

from orfaqs.lib.core.aminoacids import AminoAcid


class Protein:
    """Protein"""

    def __init__(
        self,
        amino_acid: AminoAcid | list[AminoAcid] = None,
        protein_name: str = None,
    ):
        if amino_acid is None:
            amino_acid = []
        elif isinstance(amino_acid, AminoAcid):
            amino_acid = [amino_acid]

        self._protein_name = protein_name
        self._amino_acid_chain = amino_acid
        self._sequence_complete = False

    def __str__(self) -> str:
        return ''.join(
            [amino_acid.symbol for amino_acid in self._amino_acid_chain]
        )

    @property
    def protein_name(self):
        return self._protein_name

    @property
    def amino_acid_chain(self) -> list[AminoAcid]:
        return self._amino_acid_chain

    @property
    def sequence_complete(self) -> bool:
        return self._sequence_complete

    @property
    def number_amino_acids(self) -> int:
        return len(self._amino_acid_chain)

    def add_amino_acid(self, amino_acid: AminoAcid):
        if self._sequence_complete:
            return

        self._amino_acid_chain.append(amino_acid)

    def end_sequence(self):
        self._sequence_complete = True
