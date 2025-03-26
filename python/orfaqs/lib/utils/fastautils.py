''''
FASTA Utils
'''

import os

from orfaqs.lib.core.nucleotides import (
    GenomicSequence,
    NucleotideUtils
)
from orfaqs.lib.utils.jsonutils import JsonUtils


class FASTAHeaderKeyWords:
    REF = 'ref'
    ORG = 'org'
    STRAIN = 'strain'
    MOLTYPE = 'moltype'
    CHROMOSOME = 'chromosome'
    LOCATION = 'location'
    TOP = 'top'


class FASTASequence:
    '''FASTASequence'''

    SEQUENCE_INFO_RIGHT_ARROW_DELIM = '>'
    SEQUENCE_INFO_SEMICOLON_DELIM = ';'
    _SEQUENCE_INFO_FIELD_DELIM = '['

    def __init__(self,
                 header_str: str,
                 sequence: str | list[str]):
        self._header_info = FASTASequence._parse_header(header_str)
        self._sequence_name = JsonUtils.as_json_string(
            self._header_info,
            indent=None
        )
        self._sequence = NucleotideUtils.create_sequence(
            sequence,
            self._sequence_name
        )

    @property
    def header_info(self) -> dict:
        return self._header_info

    @property
    def sequence_name(self) -> str:
        return self._sequence_name

    @property
    def sequence(self) -> GenomicSequence:
        return self._sequence

    @staticmethod
    def _parse_header(header_str: str) -> dict[str, str]:
        if not FASTAUtils.is_fasta_header(header_str):
            return None
        header_str = header_str.replace(
            FASTASequence.SEQUENCE_INFO_RIGHT_ARROW_DELIM,
            ''
        )
        header_str = header_str.replace(
            FASTASequence.SEQUENCE_INFO_SEMICOLON_DELIM,
            ''
        )

        # Grab the individual fields of the header
        header_info = {}
        header_fields = header_str.split('[')
        for field in header_fields:
            # Remove all non-essential characters
            field = field.replace('|', '')
            field = field.replace('[', ' ')
            field = field.replace(']', '')
            key = None
            if FASTAHeaderKeyWords.REF in field:
                key = FASTAHeaderKeyWords.REF
            elif FASTAHeaderKeyWords.ORG in field:
                key = FASTAHeaderKeyWords.ORG
            elif FASTAHeaderKeyWords.STRAIN in field:
                key = FASTAHeaderKeyWords.STRAIN
            elif FASTAHeaderKeyWords.MOLTYPE in field:
                key = FASTAHeaderKeyWords.MOLTYPE
            elif FASTAHeaderKeyWords.CHROMOSOME in field:
                key = FASTAHeaderKeyWords.CHROMOSOME
            elif FASTAHeaderKeyWords.LOCATION in field:
                key = FASTAHeaderKeyWords.LOCATION
            elif FASTAHeaderKeyWords.TOP in field:
                key = FASTAHeaderKeyWords.TOP

            field = field.replace(f'{key}=', '')
            field = field.replace(key, '')
            header_info[key] = field

        return header_info


class FASTAUtils:
    '''FASTAUtils'''

    _SEQUENCE_INFO_RIGHT_ARROW_DELIM = '>'
    _SEQUENCE_INFO_SEMICOLON_DELIM = ';'
    _SEQUENCE_INFO_FIELD_DELIM = '['

    @staticmethod
    def is_fasta_header(line_str: str) -> bool:
        return ((FASTAUtils._SEQUENCE_INFO_RIGHT_ARROW_DELIM in line_str) or
                (FASTAUtils._SEQUENCE_INFO_SEMICOLON_DELIM in line_str))

    @staticmethod
    def parse_file(
            file_path: str | os.PathLike) -> list[FASTASequence]:
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
                while ((line_index < len(fasta_file_lines)) and
                       not FASTAUtils.is_fasta_header(
                           fasta_file_lines[line_index])):
                    line_index += 1
                # Append a new FASTASequence object to the list
                fasta_sequences.append(
                    FASTASequence(
                        header_str,
                        fasta_file_lines[sequence_start_index:line_index]
                    )
                )

        return fasta_sequences
