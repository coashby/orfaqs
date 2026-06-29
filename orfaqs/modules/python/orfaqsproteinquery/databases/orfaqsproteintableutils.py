import enum
import logging


from orfaqs.libs.python.core.nucleotides import NucleotideUtils

_logger = logging.getLogger(__name__)


class ORFaqsProteinTableUtils:
    @staticmethod
    def uid_comment() -> str:
        return 'The unique id of the entry.'

    @staticmethod
    def source_uid_comment() -> str:
        return (
            'The identifier assigned to the sequence from '
            'which the protein was translated.'
        )

    @staticmethod
    def strand_type_comment() -> str:
        strand_types_str = NucleotideUtils.available_strand_types_str()
        return (
            f'The genomic sequence strand {strand_types_str} used during '
            'translation.'
        )

    @staticmethod
    def reading_frame_comment() -> str:
        return (
            'The reading frame (1, 2, or 3) from which the protein was '
            'transcribed.'
        )

    @staticmethod
    def orf_n_terminus_comment() -> str:
        return (
            "The base-pair position of the 5' end in the genomic sequence "
            'from which the protein was transcribed.'
        )

    @staticmethod
    def protein_comment() -> str:
        return (
            'The protein represented by its single-letter '
            'amino acid abbreviations.'
        )

    @staticmethod
    def genomic_sequence_comment() -> str:
        return (
            'The genomic sequence associated with the coding sequence using '
            'single-letter nucleotide abbreviations.'
        )

    @staticmethod
    def protein_length_comment() -> str:
        return 'The number of amino acids in the protein sequence.'


class ORFaqsProteinTableType(enum.Enum):
    DISCOVERED_PROTEINS = enum.auto()
    REFERENCE_PROTEINS = enum.auto()
