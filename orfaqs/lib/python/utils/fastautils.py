"""
FASTA Utils
"""

import enum
import logging
import os
import pathlib
import re

from orfaqs.lib.python.core.proteins import Protein
from orfaqs.lib.python.core.nucleotides import (
    DNASequence,
    RNASequence,
    NucleotideUtils,
)
from orfaqs.lib.python.core.sequence import Sequence
from orfaqs.lib.python.utils.directoryutils import DirectoryUtils
from orfaqs.lib.python.utils.jsonutils import JsonUtils

_logger = logging.getLogger(__name__)

_SUPPORTED_UID_PREFIXES: list[str] = [
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


class FASTAFileExtensionType(enum.Enum):
    FSA = enum.auto()
    FASTA = enum.auto()


class FASTASequenceType(enum.Enum):
    DNA = enum.auto()
    RNA = enum.auto()
    AMINO_ACID = enum.auto()


def _fasta_extension_type_to_str(extension_type: enum.Enum) -> str:
    return extension_type.name.lower()


def _write_fasta_str_to_file(
    file_path: pathlib.Path, fasta_str: str
) -> pathlib.Path:
    file_path = DirectoryUtils.make_path_object(file_path)
    if not FASTAUtils.is_fasta_file(file_path):
        suffix = (
            f'.{_fasta_extension_type_to_str(FASTAFileExtensionType.FASTA)}'
        )
        file_path = file_path.with_suffix(suffix)
    DirectoryUtils.mkdir_path(file_path.parent)

    with open(file_path, 'w', encoding='utf-8') as o_file:
        o_file.write(fasta_str)
        o_file.write('\n')

    return file_path


class FASTASequence:
    """FASTASequence"""

    SEQUENCE_IDENTIFIER_DELIM = '>'
    PROTEIN_SEQUENCE_STOP_CODON_DELIM = '*'
    HEADER_INFO_UID_KEY = 'uid'
    HEADER_INFO_SEQUENCE_DESCRIPTION = 'sequence_description'

    def __eq__(self, rhs):
        if not isinstance(rhs, FASTASequence):
            return False

        if (
            self.header_info != rhs.header_info
            or self.uid != rhs.uid
            or self.name != rhs.name
            or self.sequence_type != rhs.sequence_type
            or self.sequence != rhs.sequence
        ):
            return False

        return True

    def __init__(self, header_str: str, sequence_str: str):
        (self._uid, self._header_info) = FASTASequence._parse_header(
            header_str
        )
        self._name = JsonUtils.as_json_string(self._header_info, indent=None)
        sequence_str = sequence_str.replace(
            FASTASequence.PROTEIN_SEQUENCE_STOP_CODON_DELIM, ''
        )
        self._sequence = None
        try:
            self._sequence = NucleotideUtils.make_sequence_object(
                sequence=sequence_str,
                name=self._name,
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

    @staticmethod
    def _parse_header(header_str: str) -> tuple[str, dict[str, str]]:
        if not FASTAUtils.is_fasta_header(header_str):
            return None
        header_str = header_str.replace(
            FASTASequence.SEQUENCE_IDENTIFIER_DELIM, ''
        )

        # Grab the individual fields of the header
        uid = None
        header_info = {}
        header_fields = header_str.split(' ')
        number_header_fields = len(header_fields)
        if number_header_fields > 0:
            header_reference_field = header_fields[0]
            uid = re.sub(r'[<>:"/\\|?*]', '', header_reference_field)
            header_info[FASTASequence.HEADER_INFO_UID_KEY] = uid
        if number_header_fields > 1:
            header_info[FASTASequence.HEADER_INFO_SEQUENCE_DESCRIPTION] = (
                ' '.join(header_fields[1:])
            )

        return (uid, header_info)

    @property
    def header_info(self) -> dict:
        return self._header_info

    @property
    def uid(self) -> str:
        return self._uid

    @property
    def name(self) -> str:
        return self._name

    @property
    def sequence(self) -> Sequence:
        return self._sequence

    @property
    def sequence_type(self) -> FASTASequenceType:
        return self._sequence_type

    def as_fasta_str(self) -> str:
        header_info_parts: list[str] = []
        header_info_ordered_key_list = [
            FASTASequence.HEADER_INFO_UID_KEY,
            FASTASequence.HEADER_INFO_SEQUENCE_DESCRIPTION,
        ]
        for key in header_info_ordered_key_list:
            if key in self._header_info:
                header_info_parts.append(f'{self._header_info.get(key)}')

        header_info_str = ' '.join(header_info_parts)
        header_info_str = ''.join(
            [
                f'{FASTASequence.SEQUENCE_IDENTIFIER_DELIM}',
                header_info_str,
            ]
        )

        max_line_length = 60
        sequence_length = len(self._sequence)
        start_index = 0
        sequence_parts: list[str] = []
        sequence_str = self.sequence.sequence_str
        while start_index < sequence_length:
            end_index = start_index + max_line_length
            sequence_parts.append(sequence_str[start_index:end_index])
            start_index = end_index

        sequence_str = '\n'.join(sequence_parts)
        if self._sequence_type is FASTASequenceType.AMINO_ACID:
            sequence_str += FASTASequence.PROTEIN_SEQUENCE_STOP_CODON_DELIM

        return '\n'.join([header_info_str, sequence_str])

    def as_fasta_file(self, file_path: str | os.PathLike) -> pathlib.Path:
        return _write_fasta_str_to_file(
            file_path=file_path,
            fasta_str=self.as_fasta_str(),
        )


class FASTAUtils:
    """FASTAUtils"""

    @staticmethod
    def expected_file_extensions() -> list[str]:
        return [extension.name.lower() for extension in FASTAFileExtensionType]

    @staticmethod
    def is_fasta_file(file_path: str | os.PathLike) -> bool:
        file_path = DirectoryUtils.make_path_object(file_path)
        for fasta_extension in FASTAUtils.expected_file_extensions():
            suffix = f'.{fasta_extension}'
            if suffix == file_path.suffix.lower():
                return True

        return False

    @staticmethod
    def is_fasta_header(line_str: str) -> bool:
        return (len(line_str) > 0) and (
            line_str[0] == FASTASequence.SEQUENCE_IDENTIFIER_DELIM
        )

    @staticmethod
    def find_fasta_files(input_path: str | os.PathLike) -> list[pathlib.Path]:
        #######################################################################
        # Gather all FASTA files.
        input_file_paths: list[os.PathLike] = []
        if DirectoryUtils.is_file(input_path):
            input_file_paths.append(input_path)
        elif DirectoryUtils.is_directory(input_path):
            # Grab all files from the directory.
            input_file_paths = DirectoryUtils.glob_files(
                input_path, recursive=True
            )

        fasta_files: list[os.PathLike] = []
        for file_path in input_file_paths:
            if FASTAUtils.is_fasta_file(file_path):
                fasta_files.append(file_path)

        return fasta_files

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

    @staticmethod
    def export_as_fasta_file(
        file_path: str | os.PathLike,
        sequences: str | FASTASequence | Sequence | list,
        uid: str = None,
        info: str = None,
    ) -> pathlib.Path:
        input_sequences: list[any] = None
        if isinstance(sequences, list):
            input_sequences = sequences
        else:
            input_sequences = [sequences]

        fasta_str_parts: list[str] = []
        for input_sequence in input_sequences:
            if not isinstance(input_sequence, FASTASequence):
                sequence_str: str = None
                if isinstance(input_sequence, Sequence):
                    sequence_str = input_sequence.sequence_str
                    if info is None:
                        info = input_sequence.name
                elif isinstance(input_sequence, str):
                    sequence_str = input_sequence

                header_info_parts: list[str] = [
                    FASTASequence.SEQUENCE_IDENTIFIER_DELIM
                ]
                if uid is not None:
                    header_info_parts.append(uid)
                if info is not None:
                    header_info_parts.append(info)
                header_str = ' '.join(header_info_parts)
                input_sequence = FASTASequence(
                    header_str=header_str,
                    sequence_str=sequence_str,
                )

            fasta_str_parts.append(input_sequence.as_fasta_str())

        fasta_str = '\n'.join(fasta_str_parts)
        return _write_fasta_str_to_file(
            file_path=file_path,
            fasta_str=fasta_str,
        )
