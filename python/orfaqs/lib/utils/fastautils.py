"""
FASTA Utils
"""

import os
import re

from orfaqs.lib.core.nucleotides import GenomicSequence, NucleotideUtils
from orfaqs.lib.utils.jsonutils import JsonUtils

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


class FASTASequence:
    """FASTASequence"""

    ACCESSION_NUMBER_DELIM = '>'
    HEADER_INFO_ACCESSION_NUMBER_KEY = 'accession_number'
    HEADER_INFO_SEQUENCE_DESCRIPTION = 'sequence_description'

    def __init__(self, header_str: str, sequence: str | list[str]):
        (self._accession_number, self._header_info) = (
            FASTASequence._parse_header(header_str)
        )
        self._sequence_name = JsonUtils.as_json_string(
            self._header_info, indent=None
        )
        self._sequence = NucleotideUtils.create_sequence(
            sequence, self._sequence_name
        )

    @property
    def header_info(self) -> dict:
        return self._header_info

    @property
    def accession_number(self) -> str:
        return self._accession_number

    @property
    def sequence_name(self) -> str:
        return self._sequence_name

    @property
    def sequence(self) -> GenomicSequence:
        return self._sequence

    @staticmethod
    def _parse_header(header_str: str) -> tuple[str, dict[str, str]]:
        if not FASTAUtils.is_fasta_header(header_str):
            return None
        header_str = header_str.replace(
            FASTASequence.ACCESSION_NUMBER_DELIM, ''
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

    _SEQUENCE_IDENTIFIER_DELIM = '>'

    @staticmethod
    def is_fasta_header(line_str: str) -> bool:
        return (len(line_str) > 0) and (
            line_str[0] == FASTAUtils._SEQUENCE_IDENTIFIER_DELIM
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
                fasta_sequences.append(
                    FASTASequence(
                        header_str,
                        fasta_file_lines[sequence_start_index:line_index],
                    )
                )

        return fasta_sequences
