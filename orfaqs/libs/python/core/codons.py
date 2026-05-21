"""
Codons
Contains definitions for all 64 codons.
"""

import logging

from pydantic_core import core_schema
from abc import ABC, abstractmethod

from orfaqs.libs.python.core.nucleotides import (
    GenomicTriplet,
    NucleicAcid,
    NucleotideUtils,
)


_logger = logging.getLogger(__name__)


class Codon(ABC):
    """Codon"""

    _sequence: GenomicTriplet = None

    @abstractmethod
    def __init__(self, *args):
        pass

    @classmethod
    def __get_pydantic_core_schema__(cls, source_type, handler) -> str:
        return core_schema.no_info_after_validator_function(
            cls.validate, core_schema.str_schema()
        )

    def __hash__(self) -> int:
        return hash(self.sequence_str)

    def __eq__(self, rhs: any):
        if isinstance(rhs, Codon):
            return self.sequence_str == rhs.sequence_str
        elif isinstance(rhs, str):
            return self.sequence_str.upper() == rhs.upper()
        elif isinstance(rhs, list):
            if len(rhs) != Codon._NUMBER_BASES_PER_CODON:
                return False
            if isinstance(rhs[0], NucleicAcid):
                return (
                    self.sequence_str
                    == NucleotideUtils.convert_to_sequence_str(rhs)
                )

        return False

    def __str__(self) -> str:
        return self.sequence_str

    @classmethod
    def validate(cls, sequence: str) -> str:
        # No other initializer exists for Codons.
        if (
            not isinstance(sequence, str)
            or len(sequence) != GenomicTriplet.NUMBER_BASES
            or not NucleotideUtils.is_rna_sequence(sequence)
        ):
            message = '[ERROR] Invalid str input for Condon.'
            raise ValueError(message)

        return sequence

    @property
    def bases(self) -> GenomicTriplet:
        return self._sequence

    @property
    def sequence_str(self) -> str:
        return self._sequence.triplet_str

    @staticmethod
    def number_bases() -> int:
        return GenomicTriplet.NUMBER_BASES


class _UUU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UUU')


class _UUC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UUC')


class _UUA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UUA')


class _UUG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UUG')


class _UCU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UCU')


class _UCC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UCC')


class _UCA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UCA')


class _UCG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UCG')


class _UAU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UAU')


class _UAC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UAC')


class _UAA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UAA')


class _UAG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UAG')


class _UGU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UGU')


class _UGC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UGC')


class _UGA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UGA')


class _UGG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('UGG')


class _CUU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CUU')


class _CUC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CUC')


class _CUA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CUA')


class _CUG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CUG')


class _CCU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CCU')


class _CCC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CCC')


class _CCA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CCA')


class _CCG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CCG')


class _CAU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CAU')


class _CAC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CAC')


class _CAA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CAA')


class _CAG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CAG')


class _CGU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CGU')


class _CGC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CGC')


class _CGA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CGA')


class _CGG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('CGG')


class _AUU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AUU')


class _AUC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AUC')


class _AUA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AUA')


class _AUG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AUG')


class _ACU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('ACU')


class _ACC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('ACC')


class _ACA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('ACA')


class _ACG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('ACG')


class _AAU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AAU')


class _AAC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AAC')


class _AAA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AAA')


class _AAG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AAG')


class _AGU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AGU')


class _AGC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AGC')


class _AGA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AGA')


class _AGG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('AGG')


class _GUU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GUU')


class _GUC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GUC')


class _GUA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GUA')


class _GUG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GUG')


class _GCU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GCU')


class _GCC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GCC')


class _GCA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GCA')


class _GCG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GCG')


class _GAU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GAU')


class _GAC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GAC')


class _GAA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GAA')


class _GAG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GAG')


class _GGU(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GGU')


class _GGC(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GGC')


class _GGA(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GGA')


class _GGG(Codon):
    def __init__(self):
        self._sequence = GenomicTriplet('GGG')


# Constants
UUU = _UUU()
UUC = _UUC()
UUA = _UUA()
UUG = _UUG()

UCU = _UCU()
UCC = _UCC()
UCA = _UCA()
UCG = _UCG()

UAU = _UAU()
UAC = _UAC()
UAA = _UAA()
UAG = _UAG()

UGU = _UGU()
UGC = _UGC()
UGA = _UGA()
UGG = _UGG()

CUU = _CUU()
CUC = _CUC()
CUA = _CUA()
CUG = _CUG()

CCU = _CCU()
CCC = _CCC()
CCA = _CCA()
CCG = _CCG()

CAU = _CAU()
CAC = _CAC()
CAA = _CAA()
CAG = _CAG()

CGU = _CGU()
CGC = _CGC()
CGA = _CGA()
CGG = _CGG()

AUU = _AUU()
AUC = _AUC()
AUA = _AUA()
AUG = _AUG()

ACU = _ACU()
ACC = _ACC()
ACA = _ACA()
ACG = _ACG()

AAU = _AAU()
AAC = _AAC()
AAA = _AAA()
AAG = _AAG()

AGU = _AGU()
AGC = _AGC()
AGA = _AGA()
AGG = _AGG()

GUU = _GUU()
GUC = _GUC()
GUA = _GUA()
GUG = _GUG()

GCU = _GCU()
GCC = _GCC()
GCA = _GCA()
GCG = _GCG()

GAU = _GAU()
GAC = _GAC()
GAA = _GAA()
GAG = _GAG()

GGU = _GGU()
GGC = _GGC()
GGA = _GGA()
GGG = _GGG()

_AVAILABLE_CODONS: list[Codon] = [
    UUU,
    UUC,
    UUA,
    UUG,
    UCU,
    UCC,
    UCA,
    UCG,
    UAU,
    UAC,
    UAA,
    UAG,
    UGU,
    UGC,
    UGA,
    UGG,
    CUU,
    CUC,
    CUA,
    CUG,
    CCU,
    CCC,
    CCA,
    CCG,
    CAU,
    CAC,
    CAA,
    CAG,
    CGU,
    CGC,
    CGA,
    CGG,
    AUU,
    AUC,
    AUA,
    AUG,
    ACU,
    ACC,
    ACA,
    ACG,
    AAU,
    AAC,
    AAA,
    AAG,
    AGU,
    AGC,
    AGA,
    AGG,
    GUU,
    GUC,
    GUA,
    GUG,
    GCU,
    GCC,
    GCA,
    GCG,
    GAU,
    GAC,
    GAA,
    GAG,
    GGU,
    GGC,
    GGA,
    GGG,
]

_BASE_TRIPLET_CODON_LUT = {
    codon.sequence_str: codon for codon in _AVAILABLE_CODONS
}


class CodonUtils:
    @staticmethod
    def number_bases_per_codon() -> int:
        return Codon.number_bases()

    @staticmethod
    def available_codons() -> list[Codon]:
        return _AVAILABLE_CODONS

    @staticmethod
    def base_triplet_to_codon(
        base_triplet: str | GenomicTriplet,
    ):
        if isinstance(base_triplet, str):
            base_triplet = GenomicTriplet(base_triplet)

        return _BASE_TRIPLET_CODON_LUT.get(base_triplet.triplet_str)
