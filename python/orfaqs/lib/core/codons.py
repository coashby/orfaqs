'''
Codons
Contains definitions for all 64 codons.
'''
import logging

from abc import (
    ABC,
    abstractmethod
)

from orfaqs.lib.core.nucleotides import (
    NucleicAcid,
    NucleotideUtils,
    RNASequence
)


_logger = logging.getLogger(__name__)


class Codon(ABC):
    '''Codon'''
    _NUMBER_BASES_IN_CODON = 3
    _sequence: RNASequence = None

    @abstractmethod
    def __init__(self, *args):
        pass

    def __hash__(self) -> int:
        return hash(self.sequence_str)

    def __eq__(self, rhs: any):
        if isinstance(rhs, Codon):
            return self.sequence_str == rhs.sequence_str
        elif isinstance(rhs, str):
            return self.sequence_str.upper() == rhs.upper()
        elif isinstance(rhs, list):
            if (len(rhs) != Codon._NUMBER_BASES_IN_CODON):
                return False
            if isinstance(rhs[0], NucleicAcid):
                return (self.sequence_str ==
                        NucleotideUtils.convert_to_sequence_str(rhs))

        return False

    def __str__(self) -> str:
        return self.sequence_str

    @property
    def bases(self) -> RNASequence:
        return self._sequence

    @property
    def sequence_str(self) -> str:
        return self._sequence.sequence_str

    @staticmethod
    def number_bases() -> int:
        return Codon._NUMBER_BASES_IN_CODON


class _UUU(Codon):
    def __init__(self):
        self._sequence = RNASequence('UUU')


class _UUC(Codon):
    def __init__(self):
        self._sequence = RNASequence('UUC')


class _UUA(Codon):
    def __init__(self):
        self._sequence = RNASequence('UUA')


class _UUG(Codon):
    def __init__(self):
        self._sequence = RNASequence('UUG')


class _UCU(Codon):
    def __init__(self):
        self._sequence = RNASequence('UCU')


class _UCC(Codon):
    def __init__(self):
        self._sequence = RNASequence('UCC')


class _UCA(Codon):
    def __init__(self):
        self._sequence = RNASequence('UCA')


class _UCG(Codon):
    def __init__(self):
        self._sequence = RNASequence('UCG')


class _UAU(Codon):
    def __init__(self):
        self._sequence = RNASequence('UAU')


class _UAC(Codon):
    def __init__(self):
        self._sequence = RNASequence('UAC')


class _UAA(Codon):
    def __init__(self):
        self._sequence = RNASequence('UAA')


class _UAG(Codon):
    def __init__(self):
        self._sequence = RNASequence('UAG')


class _UGU(Codon):
    def __init__(self):
        self._sequence = RNASequence('UGU')


class _UGC(Codon):
    def __init__(self):
        self._sequence = RNASequence('UGC')


class _UGA(Codon):
    def __init__(self):
        self._sequence = RNASequence('UGA')


class _UGG(Codon):
    def __init__(self):
        self._sequence = RNASequence('UGG')


class _CUU(Codon):
    def __init__(self):
        self._sequence = RNASequence('CUU')


class _CUC(Codon):
    def __init__(self):
        self._sequence = RNASequence('CUC')


class _CUA(Codon):
    def __init__(self):
        self._sequence = RNASequence('CUA')


class _CUG(Codon):
    def __init__(self):
        self._sequence = RNASequence('CUG')


class _CCU(Codon):
    def __init__(self):
        self._sequence = RNASequence('CCU')


class _CCC(Codon):
    def __init__(self):
        self._sequence = RNASequence('CCC')


class _CCA(Codon):
    def __init__(self):
        self._sequence = RNASequence('CCA')


class _CCG(Codon):
    def __init__(self):
        self._sequence = RNASequence('CCG')


class _CAU(Codon):
    def __init__(self):
        self._sequence = RNASequence('CAU')


class _CAC(Codon):
    def __init__(self):
        self._sequence = RNASequence('CAC')


class _CAA(Codon):
    def __init__(self):
        self._sequence = RNASequence('CAA')


class _CAG(Codon):
    def __init__(self):
        self._sequence = RNASequence('CAG')


class _CGU(Codon):
    def __init__(self):
        self._sequence = RNASequence('CGU')


class _CGC(Codon):
    def __init__(self):
        self._sequence = RNASequence('CGC')


class _CGA(Codon):
    def __init__(self):
        self._sequence = RNASequence('CGA')


class _CGG(Codon):
    def __init__(self):
        self._sequence = RNASequence('CGG')


class _AUU(Codon):
    def __init__(self):
        self._sequence = RNASequence('AUU')


class _AUC(Codon):
    def __init__(self):
        self._sequence = RNASequence('AUC')


class _AUA(Codon):
    def __init__(self):
        self._sequence = RNASequence('AUA')


class _AUG(Codon):
    def __init__(self):
        self._sequence = RNASequence('AUG')


class _ACU(Codon):
    def __init__(self):
        self._sequence = RNASequence('ACU')


class _ACC(Codon):
    def __init__(self):
        self._sequence = RNASequence('ACC')


class _ACA(Codon):
    def __init__(self):
        self._sequence = RNASequence('ACA')


class _ACG(Codon):
    def __init__(self):
        self._sequence = RNASequence('ACG')


class _AAU(Codon):
    def __init__(self):
        self._sequence = RNASequence('AAU')


class _AAC(Codon):
    def __init__(self):
        self._sequence = RNASequence('AAC')


class _AAA(Codon):
    def __init__(self):
        self._sequence = RNASequence('AAA')


class _AAG(Codon):
    def __init__(self):
        self._sequence = RNASequence('AAG')


class _AGU(Codon):
    def __init__(self):
        self._sequence = RNASequence('AGU')


class _AGC(Codon):
    def __init__(self):
        self._sequence = RNASequence('AGC')


class _AGA(Codon):
    def __init__(self):
        self._sequence = RNASequence('AGA')


class _AGG(Codon):
    def __init__(self):
        self._sequence = RNASequence('AGG')


class _GUU(Codon):
    def __init__(self):
        self._sequence = RNASequence('GUU')


class _GUC(Codon):
    def __init__(self):
        self._sequence = RNASequence('GUC')


class _GUA(Codon):
    def __init__(self):
        self._sequence = RNASequence('GUA')


class _GUG(Codon):
    def __init__(self):
        self._sequence = RNASequence('GUG')


class _GCU(Codon):
    def __init__(self):
        self._sequence = RNASequence('GCU')


class _GCC(Codon):
    def __init__(self):
        self._sequence = RNASequence('GCC')


class _GCA(Codon):
    def __init__(self):
        self._sequence = RNASequence('GCA')


class _GCG(Codon):
    def __init__(self):
        self._sequence = RNASequence('GCG')


class _GAU(Codon):
    def __init__(self):
        self._sequence = RNASequence('GAU')


class _GAC(Codon):
    def __init__(self):
        self._sequence = RNASequence('GAC')


class _GAA(Codon):
    def __init__(self):
        self._sequence = RNASequence('GAA')


class _GAG(Codon):
    def __init__(self):
        self._sequence = RNASequence('GAG')


class _GGU(Codon):
    def __init__(self):
        self._sequence = RNASequence('GGU')


class _GGC(Codon):
    def __init__(self):
        self._sequence = RNASequence('GGC')


class _GGA(Codon):
    def __init__(self):
        self._sequence = RNASequence('GGA')


class _GGG(Codon):
    def __init__(self):
        self._sequence = RNASequence('GGG')


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
    GGG
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
    def base_triplet_to_codon(base_triplet: str | RNASequence):
        if isinstance(base_triplet, str):
            base_triplet = RNASequence(base_triplet)
        if len(base_triplet) != CodonUtils.number_bases_per_codon():
            message = ('[ERROR] Invalid base_triplet.\n'
                       f'len(\'{base_triplet}\') == {len(base_triplet)}.\n'
                       'Expected len(base_triplet) == '
                       f'{CodonUtils.number_bases_per_codon()}.\n')
            _logger.error(message)
            raise ValueError(message)

        return _BASE_TRIPLET_CODON_LUT.get(base_triplet.sequence_str)
