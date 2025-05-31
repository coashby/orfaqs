"""
Abstraction of the 22 naturally occurring amino acids.
"""

from abc import ABC, abstractmethod


class AminoAcid(ABC):
    """AminoAcid"""

    def __hash__(self):
        return hash(self.abbreviation)

    def __eq__(self, rhs: any) -> bool:
        if isinstance(rhs, AminoAcid):
            return rhs.abbreviation == self.abbreviation

        return False

    def __str__(self) -> str:
        return self.symbol

    @property
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    @abstractmethod
    def abbreviation(self) -> str:
        pass

    @property
    @abstractmethod
    def symbol(self) -> str:
        pass


class _Alanine(AminoAcid):
    """_Alanine"""

    @property
    def name(self) -> str:
        return 'alanine'

    @property
    def abbreviation(self) -> str:
        return 'Ala'

    @property
    def symbol(self) -> str:
        return 'A'


class _Arginine(AminoAcid):
    """_Arginine"""

    @property
    def name(self) -> str:
        return 'arginine'

    @property
    def abbreviation(self) -> str:
        return 'Arg'

    @property
    def symbol(self) -> str:
        return 'R'


class _Asparagine(AminoAcid):
    """_Asparagine"""

    @property
    def name(self) -> str:
        return 'asparagine'

    @property
    def abbreviation(self) -> str:
        return 'Asn'

    @property
    def symbol(self) -> str:
        return 'N'


class _Aspartate(AminoAcid):
    """_Aspartate"""

    @property
    def name(self) -> str:
        return 'aspartate'

    @property
    def abbreviation(self) -> str:
        return 'Asp'

    @property
    def symbol(self) -> str:
        return 'D'


class _Cysteine(AminoAcid):
    """_Cysteine"""

    @property
    def name(self) -> str:
        return 'cysteine'

    @property
    def abbreviation(self) -> str:
        return 'Cys'

    @property
    def symbol(self) -> str:
        return 'C'


class _Glutamine(AminoAcid):
    """_Glutamine"""

    @property
    def name(self) -> str:
        return 'glutamine'

    @property
    def abbreviation(self) -> str:
        return 'Gln'

    @property
    def symbol(self) -> str:
        return 'Q'


class _Glutamate(AminoAcid):
    """_Glutamate"""

    @property
    def name(self) -> str:
        return 'glutamate'

    @property
    def abbreviation(self) -> str:
        return 'Glu'

    @property
    def symbol(self) -> str:
        return 'E'


class _Glycine(AminoAcid):
    """_Glycine"""

    @property
    def name(self) -> str:
        return 'glycine'

    @property
    def abbreviation(self) -> str:
        return 'Gly'

    @property
    def symbol(self) -> str:
        return 'G'


class _Histidine(AminoAcid):
    """_Histidine"""

    @property
    def name(self) -> str:
        return 'histidine'

    @property
    def abbreviation(self) -> str:
        return 'His'

    @property
    def symbol(self) -> str:
        return 'H'


class _Isoleucine(AminoAcid):
    """_Isoleucine"""

    @property
    def name(self) -> str:
        return 'isoleucine'

    @property
    def abbreviation(self) -> str:
        return 'Ile'

    @property
    def symbol(self) -> str:
        return 'I'


class _Leucine(AminoAcid):
    """_Leucine"""

    @property
    def name(self) -> str:
        return 'leucine'

    @property
    def abbreviation(self) -> str:
        return 'Leu'

    @property
    def symbol(self) -> str:
        return 'L'


class _Lysine(AminoAcid):
    """_Lysine"""

    @property
    def name(self) -> str:
        return 'lysine'

    @property
    def abbreviation(self) -> str:
        return 'Lys'

    @property
    def symbol(self) -> str:
        return 'K'


class _Methionine(AminoAcid):
    """_Methionine"""

    @property
    def name(self) -> str:
        return 'methionine'

    @property
    def abbreviation(self) -> str:
        return 'Met'

    @property
    def symbol(self) -> str:
        return 'M'


class _Phenylalanine(AminoAcid):
    """_Phenylalanine"""

    @property
    def name(self) -> str:
        return 'phenylalanine'

    @property
    def abbreviation(self) -> str:
        return 'Phe'

    @property
    def symbol(self) -> str:
        return 'F'


class _Proline(AminoAcid):
    """_Proline"""

    @property
    def name(self) -> str:
        return 'proline'

    @property
    def abbreviation(self) -> str:
        return 'Pro'

    @property
    def symbol(self) -> str:
        return 'P'


class _Serine(AminoAcid):
    """_Serine"""

    @property
    def name(self) -> str:
        return 'serine'

    @property
    def abbreviation(self) -> str:
        return 'Ser'

    @property
    def symbol(self) -> str:
        return 'S'


class _Threonine(AminoAcid):
    """_Threonine"""

    @property
    def name(self) -> str:
        return 'threonine'

    @property
    def abbreviation(self) -> str:
        return 'Thr'

    @property
    def symbol(self) -> str:
        return 'T'


class _Tryptophan(AminoAcid):
    """_Tryptophan"""

    @property
    def name(self) -> str:
        return 'tryptophan'

    @property
    def abbreviation(self) -> str:
        return 'Trp'

    @property
    def symbol(self) -> str:
        return 'W'


class _Tyrosine(AminoAcid):
    """_Tyrosine"""

    @property
    def name(self) -> str:
        return 'tyrosine'

    @property
    def abbreviation(self) -> str:
        return 'Tyr'

    @property
    def symbol(self) -> str:
        return 'Y'


class _Valine(AminoAcid):
    """_Valine"""

    @property
    def name(self) -> str:
        return 'valine'

    @property
    def abbreviation(self) -> str:
        return 'Val'

    @property
    def symbol(self) -> str:
        return 'V'


class _Selenocysteine(AminoAcid):
    """_Selenocysteine"""

    @property
    def name(self) -> str:
        return 'selenocysteine'

    @property
    def abbreviation(self) -> str:
        return 'Sec'

    @property
    def symbol(self) -> str:
        return 'U'


class _Pyrrolysine(AminoAcid):
    """_Pyrrolysine"""

    @property
    def name(self) -> str:
        return 'pyrrolysine'

    @property
    def abbreviation(self) -> str:
        return 'Pyl'

    @property
    def symbol(self) -> str:
        return 'O'


# Constants
ALANINE = _Alanine()
ARGININE = _Arginine()
ASPARAGINE = _Asparagine()
ASPARTATE = _Aspartate()
CYSTEINE = _Cysteine()
GLUTAMINE = _Glutamine()
GLUTAMATE = _Glutamate()
GLYCINE = _Glycine()
HISTIDINE = _Histidine()
ISOLEUCINE = _Isoleucine()
LEUCINE = _Leucine()
LYSINE = _Lysine()
METHIONINE = _Methionine()
PHENYLALANINE = _Phenylalanine()
PROLINE = _Proline()
SERINE = _Serine()
THREONINE = _Threonine()
TRYPTOPHAN = _Tryptophan()
TYROSINE = _Tyrosine()
VALINE = _Valine()
SELENOCYSTEINE = _Selenocysteine()
PYRROLYSINE = _Pyrrolysine()

_AVAILABLE_AMINO_ACIDS: list[AminoAcid] = [
    ALANINE,
    ARGININE,
    ASPARAGINE,
    ASPARTATE,
    CYSTEINE,
    GLUTAMINE,
    GLUTAMATE,
    GLYCINE,
    HISTIDINE,
    ISOLEUCINE,
    LEUCINE,
    LYSINE,
    METHIONINE,
    PHENYLALANINE,
    PROLINE,
    SERINE,
    THREONINE,
    TRYPTOPHAN,
    TYROSINE,
    VALINE,
    SELENOCYSTEINE,
    PYRROLYSINE,
]

_AMINO_ACID_SYMBOL_LUT = {
    amino_acid.symbol: amino_acid for amino_acid in _AVAILABLE_AMINO_ACIDS
}


class AminoAcidUtils:
    """AminoAcidUtils"""

    @staticmethod
    def available_amino_acids() -> list[AminoAcid]:
        return _AVAILABLE_AMINO_ACIDS

    @staticmethod
    def symbol_to_amino_acid(symbol: str) -> AminoAcid:
        return _AMINO_ACID_SYMBOL_LUT.get(symbol)
