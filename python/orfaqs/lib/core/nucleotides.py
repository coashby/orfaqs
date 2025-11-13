"""
Nucleotides
"""

import enum
import logging
from abc import ABC, abstractmethod

_logger = logging.getLogger(__name__)


class NucleicAcid(ABC):
    """NucleicAcid"""

    def __init__(self):
        self._convert_to_str_using_name = False

    def __hash__(self) -> int:
        return hash(self.name)

    def __eq__(self, rhs: any) -> bool:
        if isinstance(rhs, str):
            return (self.name == rhs.lower()) or self.symbol == rhs.upper()
        elif isinstance(rhs, NucleicAcid):
            return (self.name == rhs.name) and (self.symbol == rhs.symbol)

        return False

    def __str__(self) -> str:
        if self._convert_to_str_using_name:
            return self.name
        return self.symbol

    @property
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    @abstractmethod
    def symbol(self) -> str:
        pass

    @property
    def convert_to_str_using_name(self) -> bool:
        return self._convert_to_str_using_name

    @convert_to_str_using_name.setter
    def convert_to_str_using_name(self, value: bool):
        self._convert_to_str_using_name = value


class _Adenine(NucleicAcid):
    """_Adenine"""

    @property
    def name(self) -> str:
        return 'adenine'

    @property
    def symbol(self) -> str:
        return 'A'


class _Cytosine(NucleicAcid):
    """_Cytosine"""

    @property
    def name(self) -> str:
        return 'cytosine'

    @property
    def symbol(self) -> str:
        return 'C'


class _Guanine(NucleicAcid):
    """_Guanine"""

    @property
    def name(self) -> str:
        return 'guanine'

    @property
    def symbol(self) -> str:
        return 'G'


class _Thymine(NucleicAcid):
    """_Thymine"""

    @property
    def name(self) -> str:
        return 'thymine'

    @property
    def symbol(self) -> str:
        return 'T'


class _Uracil(NucleicAcid):
    """_Uracil"""

    @property
    def name(self) -> str:
        return 'uracil'

    @property
    def symbol(self) -> str:
        return 'U'


ADENINE = _Adenine()
CYTOSINE = _Cytosine()
GUANINE = _Guanine()
THYMINE = _Thymine()
URACIL = _Uracil()

_NUCLEIC_ACID_SYMBOL_LUT = {
    ADENINE.symbol: ADENINE,
    CYTOSINE.symbol: CYTOSINE,
    GUANINE.symbol: GUANINE,
    THYMINE.symbol: THYMINE,
    URACIL.symbol: URACIL,
}

_NUCLEIC_ACID_NAME_LUT = {
    ADENINE.name: ADENINE,
    CYTOSINE.name: CYTOSINE,
    GUANINE.name: GUANINE,
    THYMINE.name: THYMINE,
    URACIL.name: URACIL,
}


class StrandType(enum.Enum):
    NEGATIVE_STRAND = 'negative_strand'
    POSITIVE_STRAND = 'positive_strand'
    SINGLE_STRANDED = 'single_stranded'

    def __str__(self) -> str:
        return self.value


class GenomicSequence:
    """GenomicSequence"""

    def __init__(
        self,
        sequence: (str | list[str] | list[NucleicAcid]),
        strand_type: StrandType = None,
        name: str = None,
    ):
        self._sequence_str: str
        self._name = name
        self._strand_type = strand_type

        if isinstance(sequence, GenomicSequence):
            self._sequence_str = sequence._sequence_str
            if self._name is None:
                self._name = sequence._name
            if self._strand_type is None:
                self._strand_type = sequence._strand_type
        else:
            if isinstance(sequence, list):
                sequence = ''.join(sequence)

            sequence = sequence.lower()
            for base in sequence:
                if base not in list(self.available_bases()):
                    message = (
                        f'[ERROR] The following symbol {base} is not '
                        "allowed in GenomicSequence object's of type: "
                        f'{self.__class__.__name__}'
                    )
                    _logger.error(message)
                    raise ValueError(message)

            self._sequence_str = sequence

        if self._strand_type is None:
            self._strand_type = StrandType.POSITIVE_STRAND

        self._base_iterator = self._base_generator()

    def __getitem__(self, index) -> any:
        item = self._sequence_str[index]
        if len(item) > 1:
            return self.__class__(item, self._strand_type, self._name)

        return _NUCLEIC_ACID_SYMBOL_LUT[item]

    def __contains__(self, region: any) -> bool:
        if isinstance(region, str):
            return region in self._sequence_str
        elif isinstance(region, GenomicSequence):
            return region._sequence_str in self._sequence_str
        elif isinstance(region, list):
            return self.__class__(region)._sequence_str in self._sequence_str

        return False

    def __str__(self) -> str:
        return self._sequence_str

    def __len__(self) -> int:
        return len(self._sequence_str)

    @staticmethod
    @abstractmethod
    def available_bases() -> list[NucleicAcid]:
        pass

    @staticmethod
    @abstractmethod
    def base_complement(base: str | NucleicAcid) -> NucleicAcid:
        pass

    @property
    def name(self) -> str:
        return self._name

    @property
    def strand_type(self) -> str:
        return self._strand_type

    @property
    def sequence_str(self) -> str:
        return self._sequence_str

    @property
    def sequence_length(self) -> int:
        return len(self._sequence_str)

    @property
    def base_iterator(self):
        return self._base_generator

    def _base_generator(self):
        for base in self._sequence_str:
            yield base

    def reset_base_iterator(self):
        self._base_iterator = self._base_generator()


class DNASequence(GenomicSequence):
    """DNASequence"""

    @staticmethod
    def available_bases() -> list[NucleicAcid]:
        return list([ADENINE, THYMINE, CYTOSINE, GUANINE])

    @staticmethod
    def base_complement(base: str | NucleicAcid) -> NucleicAcid:
        if base not in DNASequence.available_bases():
            message = (
                f'[ERROR] {base} is not recognized as a nucleic acid for DNA'
            )
            _logger.error(message)
            raise ValueError(message)

        if base == ADENINE:
            return THYMINE
        elif base == THYMINE:
            return ADENINE
        elif base == CYTOSINE:
            return GUANINE
        elif base == GUANINE:
            return CYTOSINE

        return None


class RNASequence(GenomicSequence):
    """RNASequence"""

    @staticmethod
    def available_bases() -> list[NucleicAcid]:
        return [ADENINE, URACIL, CYTOSINE, GUANINE]

    @staticmethod
    def base_complement(base: str | NucleicAcid) -> NucleicAcid:
        if base not in RNASequence.available_bases():
            message = (
                f'[ERROR] {base} is not recognized as a nucleic acid for RNA'
            )
            _logger.error(message)
            raise ValueError(message)

        if base == ADENINE:
            return URACIL
        elif base == URACIL:
            return ADENINE
        elif base == CYTOSINE:
            return GUANINE
        elif base == GUANINE:
            return CYTOSINE

        return None


class NucleotideUtils:
    """NucleotideUtils"""

    @staticmethod
    def available_strand_types() -> list[StrandType]:
        return [strand_type for strand_type in StrandType]

    @staticmethod
    def available_strand_types_str() -> list[str]:
        return [strand_type.value for strand_type in StrandType]

    @staticmethod
    def create_sequence(
        sequence: (str | list[str] | list[NucleicAcid] | GenomicSequence),
        strand_type: StrandType = None,
        name: str = None,
    ) -> GenomicSequence:
        try:
            return DNASequence(sequence, strand_type, name)
        except ValueError:
            try:
                return RNASequence(sequence, strand_type, name)
            except ValueError as e:
                message = (
                    '[ERROR] The sequence could not be interpreted as '
                    'a DNA sequence or an RNA sequence. Check the '
                    'input for errors.'
                )
                _logger.error(message)
                raise ValueError(message) from e

    @staticmethod
    def is_dna_sequence(
        sequence: (str | list[str] | list[NucleicAcid]),
    ) -> bool:
        if isinstance(sequence, list):
            sequence = ''.join(sequence)

        for base in sequence:
            if base not in list(DNASequence.available_bases()):
                return False

        return True

    @staticmethod
    def is_rna_sequence(
        sequence: (str | list[str] | list[NucleicAcid]),
    ) -> bool:
        if isinstance(sequence, list):
            sequence = ''.join(sequence)

        for base in sequence:
            if base not in list(RNASequence.available_bases()):
                return False

        return True

    @staticmethod
    def dna_nucleotides() -> list[NucleicAcid]:
        return list([ADENINE, THYMINE, GUANINE, CYTOSINE])

    @staticmethod
    def rna_nucleotides() -> list[NucleicAcid]:
        return list([ADENINE, URACIL, GUANINE, CYTOSINE])

    @staticmethod
    def convert_to_rna_base(base: NucleicAcid) -> NucleicAcid:
        if base == THYMINE:
            return URACIL
        return base

    @staticmethod
    def convert_to_dna_base(base: NucleicAcid) -> NucleicAcid:
        if base == URACIL:
            return THYMINE

        return base

    @staticmethod
    def symbol_to_nucleic_acid(symbol: str) -> NucleicAcid:
        return _NUCLEIC_ACID_SYMBOL_LUT.get(symbol)

    @staticmethod
    def name_to_nucleic_acid(name: str) -> NucleicAcid:
        return _NUCLEIC_ACID_NAME_LUT.get(name)

    @staticmethod
    def convert_to_sequence_str(bases: list[NucleicAcid]) -> str:
        symbol_list: list[str] = []
        for base in bases:
            symbol_list.append(base.symbol)

        return ''.join(symbol_list)

    @staticmethod
    def reverse_str(genomic_sequence: GenomicSequence) -> str:
        return genomic_sequence.sequence_str[::-1]

    @staticmethod
    def complement_str(genomic_sequence: GenomicSequence) -> str:
        complement_str = ''
        for base in genomic_sequence.sequence_str:
            complement_str += genomic_sequence.base_complement(base).symbol

        return complement_str

    @staticmethod
    def reverse_complement_str(genomic_sequence: GenomicSequence) -> str:
        return NucleotideUtils.complement_str(genomic_sequence)[::-1]

    @staticmethod
    def reverse_complement(
        genomic_sequence: GenomicSequence,
    ) -> GenomicSequence:
        reverse_complement_sequence = NucleotideUtils.reverse_complement_str(
            genomic_sequence
        )
        name = genomic_sequence.name
        strand_type = genomic_sequence.strand_type
        if strand_type is StrandType.POSITIVE_STRAND:
            strand_type = StrandType.NEGATIVE_STRAND
        elif strand_type is StrandType.NEGATIVE_STRAND:
            strand_type = StrandType.POSITIVE_STRAND

        if isinstance(genomic_sequence, DNASequence):
            return DNASequence(
                sequence=reverse_complement_sequence,
                name=name,
                strand_type=strand_type,
            )
        elif isinstance(genomic_sequence, RNASequence):
            return RNASequence(
                sequence=reverse_complement_sequence,
                name=name,
                strand_type=strand_type,
            )
        return None

    @staticmethod
    def find_sequence(
        genomic_sequence: GenomicSequence, sequence_pattern: any
    ) -> int:
        sequence_pattern_str = None
        if isinstance(sequence_pattern, str):
            sequence_pattern_str = sequence_pattern
        elif isinstance(sequence_pattern, GenomicSequence):
            sequence_pattern_str = sequence_pattern.sequence_str
        elif isinstance(sequence_pattern, list):
            sequence_pattern_str = genomic_sequence.__class__(
                sequence_pattern
            ).sequence_str

        if isinstance(sequence_pattern_str, str):
            return genomic_sequence.sequence_str(sequence_pattern_str)

        return None
