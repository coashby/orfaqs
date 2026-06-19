"""
Nucleotides
"""

import enum
import logging

from abc import ABC, abstractmethod

from orfaqs.libs.python.core.sequence import Sequence

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

_NUCLEOTIDES = [
    ADENINE,
    CYTOSINE,
    GUANINE,
    THYMINE,
    URACIL,
]
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


class StrandType(enum.StrEnum):
    NEGATIVE_STRAND = 'negative_strand'
    POSITIVE_STRAND = 'positive_strand'
    SINGLE_STRANDED = 'single_stranded'

    def __str__(self) -> str:
        return self.value


class GenomicTriplet:
    NUMBER_BASES = 3

    def __eq__(self, value):
        value = str(value)
        return self._triplet_str == value.upper()

    def __init__(self, *nucleotides):
        raise_error = False
        nucleotide_str = ''
        if len(nucleotides) == 1:
            nucleotides = nucleotides[0]
            if isinstance(nucleotides, str) or isinstance(
                nucleotides, GenomicSequence
            ):
                nucleotide_str = str(nucleotides)
            else:
                raise_error = True

        elif len(nucleotides) == GenomicTriplet.NUMBER_BASES:
            if all(
                isinstance(nucleotide, NucleicAcid)
                for nucleotide in nucleotides
            ) or all(
                isinstance(nucleotide, str) for nucleotide in nucleotides
            ):
                for nucleotide in nucleotides:
                    nucleotide_str += str(nucleotide)
            else:
                raise_error = True

        if len(nucleotide_str) > GenomicTriplet.NUMBER_BASES:
            raise_error = True

        nucleotide_str = nucleotide_str.upper()
        nucleotide_symbols = [str(nucleotide) for nucleotide in _NUCLEOTIDES]
        for nucleotide in nucleotide_str:
            if nucleotide not in nucleotide_symbols:
                raise_error = True
                break

        if str(URACIL) in nucleotide_str and str(THYMINE) in nucleotide_str:
            raise_error = True

        if raise_error:
            message = '[ERROR] Unsupported input type.'
            _logger.error(message)
            raise ValueError(message)

        self._triplet_str = nucleotide_str.upper()

    @property
    def triplet_str(self) -> str:
        return self._triplet_str


class GenomicSequence(Sequence):
    """GenomicSequence"""

    def __init__(
        self,
        sequence: (str | list[str] | list[NucleicAcid]),
        strand_type: StrandType = None,
        name: str = None,
        uid: str = None,
        info: str = None,
        log_errors: bool = True,
        raise_errors: bool = True,
    ):
        super().__init__(
            sequence=sequence,
            name=name,
            uid=uid,
            info=info,
            log_errors=log_errors,
            raise_errors=raise_errors,
        )
        self._strand_type = strand_type

        if isinstance(sequence, GenomicSequence):
            self._sequence_str = sequence._sequence_str
            if self._name is None:
                self._name = sequence._name
            if self._strand_type is None:
                self._strand_type = sequence._strand_type

        if self._strand_type is None:
            self._strand_type = StrandType.POSITIVE_STRAND

    def __getitem__(self, index) -> any:
        item = self._sequence_str[index]
        if len(item) > 1:
            return self.__class__(item, self._strand_type, self._name)

        return _NUCLEIC_ACID_SYMBOL_LUT[item]

    def __len__(self) -> int:
        return len(self._sequence_str)

    @classmethod
    def available_symbols(cls) -> list[str]:
        return [base.symbol for base in cls.available_bases()]

    @staticmethod
    @abstractmethod
    def available_bases() -> list[NucleicAcid]:
        pass

    @staticmethod
    @abstractmethod
    def base_complement(base: str | NucleicAcid) -> NucleicAcid:
        pass

    @staticmethod
    def number_bases_per_triplet() -> int:
        return GenomicTriplet.NUMBER_BASES

    def number_triplets(self) -> int:
        return int(len(self._sequence_str) / GenomicTriplet.NUMBER_BASES)

    @property
    def strand_type(self) -> StrandType:
        return self._strand_type

    def triplet(self, index) -> GenomicTriplet | None:
        base_index = index * GenomicTriplet.NUMBER_BASES
        end_base_index = base_index + GenomicTriplet.NUMBER_BASES
        if end_base_index > len(self._sequence_str):
            return None
        triplet = self._sequence_str[base_index:end_base_index]
        return GenomicTriplet(triplet)


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
    def make_sequence_object(
        sequence: (str | list[str] | list[NucleicAcid] | GenomicSequence),
        strand_type: StrandType = None,
        name: str = None,
        raise_errors: bool = True,
    ) -> GenomicSequence:
        try:
            return DNASequence(
                sequence,
                strand_type,
                name,
                log_errors=False,
                raise_errors=raise_errors,
            )
        except ValueError:
            try:
                return RNASequence(
                    sequence,
                    strand_type,
                    name,
                    log_errors=False,
                    raise_errors=raise_errors,
                )
            except ValueError as e:
                if raise_errors:
                    raise ValueError from e

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
