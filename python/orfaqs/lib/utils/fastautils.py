"""
FASTA Utils
"""

import enum
import logging
import os
import re

from orfaqs.lib.core.proteins import Protein
from orfaqs.lib.core.nucleotides import (
    DNASequence,
    RNASequence,
    NucleotideUtils,
)
from orfaqs.lib.core.sequence import Sequence
from orfaqs.lib.utils.jsonutils import JsonUtils

_logger = logging.getLogger(__name__)

_SUPPORTED_ACCESSION_NUMBER_PREFIXES: list[str] = [
    'NC_',
    'NG_',
    'NM_',
    'NP_',
    'NR_',
    'NT_',
    'NW_',
    'XM_',
    'XP_',
    'XR_',
]


class FASTASequenceType(enum.Enum):
    DNA = enum.auto()
    RNA = enum.auto()
    AMINO_ACID = enum.auto()


class FASTASequence:
    """FASTASequence"""

    SEQUENCE_IDENTIFIER_DELIM = '>'
    PROTEIN_SEQUENCE_STOP_CODON_DELIM = '*'
    HEADER_INFO_ACCESSION_NUMBER_KEY = 'accession_number'
    HEADER_INFO_SEQUENCE_DESCRIPTION = 'sequence_description'

    def __init__(self, header_str: str, sequence_str: str):
        (self._accession_number, self._header_info) = (
            FASTASequence._parse_header(header_str)
        )
        self._name = JsonUtils.as_json_string(self._header_info, indent=None)
        sequence_str = sequence_str.replace(
            FASTASequence.PROTEIN_SEQUENCE_STOP_CODON_DELIM, ''
        )
        self._sequence = None
        try:
            self._sequence = NucleotideUtils.create_sequence(
                sequence=sequence_str, name=self._name
            )
        except ValueError:
            try:
                self._sequence = Protein(
                    sequence=sequence_str,
                    name=self.name,
                    sequence_complete=True,
                )
            except ValueError as e:
                message = (
                    '[ERROR] No Sequence type supports the provided '
                    'sequence string. Supported types are: genomic '
                    'sequence(RNA and DNA) and protein sequences.\n'
                    f'{e}'
                )
                _logger.error(message)
                raise ValueError(message) from e

        self._sequence_type: FASTASequenceType = None
        if isinstance(self._sequence, DNASequence):
            self._sequence_type = FASTASequenceType.DNA
        elif isinstance(self._sequence, RNASequence):
            self._sequence_type = FASTASequenceType.RNA
        elif isinstance(self._sequence, Protein):
            self._sequence_type = FASTASequenceType.AMINO_ACID

    @property
    def header_info(self) -> dict:
        return self._header_info

    @property
    def accession_number(self) -> str:
        return self._accession_number

    @property
    def name(self) -> str:
        return self._name

    @property
    def sequence(self) -> Sequence:
        return self._sequence

    @property
    def sequence_type(self) -> FASTASequenceType:
        return self._sequence_type

    @staticmethod
    def _parse_header(header_str: str) -> tuple[str, dict[str, str]]:
        if not FASTAUtils.is_fasta_header(header_str):
            return None
        header_str = header_str.replace(
            FASTASequence.SEQUENCE_IDENTIFIER_DELIM, ''
        )

        # Grab the individual fields of the header
        accession_number = None
        header_info = {}
        header_fields = header_str.split(' ')
        number_header_fields = len(header_fields)
        if number_header_fields > 0:
            header_reference_field = header_fields[0]
            for prefix in _SUPPORTED_ACCESSION_NUMBER_PREFIXES:
                if prefix in header_reference_field:
                    accession_number_fields = header_reference_field.split(
                        prefix
                    )
                    # Remove all non-alphanumeric characters from the
                    # accession_number.
                    accession_number = re.sub(
                        r'[^a-zA-Z0-9]', '', accession_number_fields[-1]
                    )
                    accession_number = f'{prefix}{accession_number}'

            header_info[FASTASequence.HEADER_INFO_ACCESSION_NUMBER_KEY] = (
                accession_number
            )
        if number_header_fields > 1:
            header_info[FASTASequence.HEADER_INFO_SEQUENCE_DESCRIPTION] = (
                ' '.join(header_fields[1:])
            )

        return (accession_number, header_info)


class FASTAUtils:
    """FASTAUtils"""

    @staticmethod
    def is_fasta_header(line_str: str) -> bool:
        return (len(line_str) > 0) and (
            line_str[0] == FASTASequence.SEQUENCE_IDENTIFIER_DELIM
        )

    @staticmethod
    def parse_file(file_path: str | os.PathLike) -> list[FASTASequence]:
        fasta_file_lines = []
        with open(file_path, 'r', encoding='utf-8') as i_file:
            for line in i_file:
                fasta_file_lines.append(line.strip())

        fasta_sequences: list[FASTASequence] = []
        line_index = 0
        while line_index < len(fasta_file_lines):
            if FASTAUtils.is_fasta_header(fasta_file_lines[line_index]):
                header_str = fasta_file_lines[line_index]
                line_index += 1
                sequence_start_index = line_index
                while (
                    line_index < len(fasta_file_lines)
                ) and not FASTAUtils.is_fasta_header(
                    fasta_file_lines[line_index]
                ):
                    line_index += 1
                # Append a new FASTASequence object to the list
                sequence_str = ''.join(
                    fasta_file_lines[sequence_start_index:line_index]
                )
                fasta_sequences.append(
                    FASTASequence(
                        header_str,
                        sequence_str,
                    )
                )

        return fasta_sequences
