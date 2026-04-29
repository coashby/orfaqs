"""
Fasta Utils test module.
"""

import pathlib
import pytest

from pytest_datadir.plugin import LazyDataDir

from orfaqs.lib.python.utils.fastautils import (
    FASTASequenceType,
    FASTASequence,
    FASTAUtils,
)


@pytest.fixture
def inputs_directory_path(lazy_shared_datadir: LazyDataDir) -> pathlib.Path:
    return lazy_shared_datadir.original_datadir / 'inputs'


@pytest.fixture
def dna_sequences_fasta_file_path(
    inputs_directory_path: pathlib.Path,
) -> pathlib.Path:
    return inputs_directory_path / 'dna-sequences.fasta'


@pytest.fixture
def rna_sequences_fasta_file_path(
    inputs_directory_path: pathlib.Path,
) -> pathlib.Path:
    return inputs_directory_path / 'rna-sequences.fasta'


@pytest.fixture
def protein_sequences_fasta_file_path(
    inputs_directory_path: pathlib.Path,
) -> pathlib.Path:
    return inputs_directory_path / 'protein-sequences.fasta'


@pytest.fixture
def invalid_sequence_fasta_file_path(
    inputs_directory_path: pathlib.Path,
) -> pathlib.Path:
    return inputs_directory_path / 'invalid-sequence.fasta'


class TestFastaUtils:
    @staticmethod
    def _validate_number_sequences(
        file_path: str, fasta_sequences: list[FASTASequence]
    ):
        expected_number_sequences = None
        with open(file_path, mode='r', encoding='utf-8') as i_file:
            expected_number_sequences = i_file.read().count(
                FASTASequence.SEQUENCE_IDENTIFIER_DELIM
            )
        assert len(fasta_sequences) == expected_number_sequences

    @staticmethod
    def test_read_dna_sequences(dna_sequences_fasta_file_path: pathlib.Path):
        fasta_sequences = FASTAUtils.parse_file(dna_sequences_fasta_file_path)
        TestFastaUtils._validate_number_sequences(
            dna_sequences_fasta_file_path,
            fasta_sequences,
        )
        for fasta_sequence in fasta_sequences:
            assert fasta_sequence.sequence_type == FASTASequenceType.DNA

    @staticmethod
    def test_read_rna_sequences(rna_sequences_fasta_file_path: pathlib.Path):
        fasta_sequences = FASTAUtils.parse_file(rna_sequences_fasta_file_path)
        TestFastaUtils._validate_number_sequences(
            rna_sequences_fasta_file_path,
            fasta_sequences,
        )
        for fasta_sequence in fasta_sequences:
            assert fasta_sequence.sequence_type == FASTASequenceType.RNA

    @staticmethod
    def test_read_protein_sequences(
        protein_sequences_fasta_file_path: pathlib.Path,
    ):
        fasta_sequences = FASTAUtils.parse_file(
            protein_sequences_fasta_file_path
        )
        TestFastaUtils._validate_number_sequences(
            protein_sequences_fasta_file_path,
            fasta_sequences,
        )
        for fasta_sequence in fasta_sequences:
            assert fasta_sequence.sequence_type == FASTASequenceType.AMINO_ACID

    @staticmethod
    def test_invalid_sequence(invalid_sequence_fasta_file_path: pathlib.Path):
        value_error_exception_thrown = False
        try:
            _ = FASTAUtils.parse_file(invalid_sequence_fasta_file_path)
        except ValueError:
            value_error_exception_thrown = True

        assert value_error_exception_thrown

    @staticmethod
    def _validate_single_sequence_exported_fasta_file(
        exported_fasta_file_path: pathlib.Path,
        expected_fasta_sequence: FASTASequence,
    ):
        fasta_sequences = FASTAUtils.parse_file(exported_fasta_file_path)
        assert len(fasta_sequences) == 1
        fasta_sequence = fasta_sequences[0]

        assert fasta_sequence == expected_fasta_sequence

    @staticmethod
    def _validate_export_single_sequence_fasta_file_from_fasta_sequences(
        fasta_file_path: pathlib.Path,
        export_directory_path: pathlib.Path,
    ):
        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        for index, fasta_sequence in enumerate(fasta_sequences):
            file_name = f'exported-sequence-{index}'
            file_path = export_directory_path / file_name
            file_path = fasta_sequence.as_fasta_file(file_path)
            # Validate the exported sequence
            TestFastaUtils._validate_single_sequence_exported_fasta_file(
                exported_fasta_file_path=file_path,
                expected_fasta_sequence=fasta_sequence,
            )

    @staticmethod
    def test_export_single_sequence_fasta_file_from_fasta_dna_sequences(
        dna_sequences_fasta_file_path: pathlib.Path,
        tmp_path: pathlib.Path,
    ):
        TestFastaUtils._validate_export_single_sequence_fasta_file_from_fasta_sequences(
            fasta_file_path=dna_sequences_fasta_file_path,
            export_directory_path=tmp_path,
        )

    @staticmethod
    def test_export_single_sequence_fasta_file_from_fasta_rna_sequence(
        rna_sequences_fasta_file_path: pathlib.Path,
        tmp_path: pathlib.Path,
    ):
        TestFastaUtils._validate_export_single_sequence_fasta_file_from_fasta_sequences(
            fasta_file_path=rna_sequences_fasta_file_path,
            export_directory_path=tmp_path,
        )

    @staticmethod
    def test_export_single_sequence_fasta_file_from_fasta_protein_sequences(
        protein_sequences_fasta_file_path: pathlib.Path,
        tmp_path: pathlib.Path,
    ):
        TestFastaUtils._validate_export_single_sequence_fasta_file_from_fasta_sequences(
            fasta_file_path=protein_sequences_fasta_file_path,
            export_directory_path=tmp_path,
        )

    @staticmethod
    def _validate_multi_sequence_exported_fasta_file(
        exported_fasta_file_path: pathlib.Path,
        expected_fasta_file_path: pathlib.Path,
    ):
        fasta_sequences = FASTAUtils.parse_file(exported_fasta_file_path)
        expected_fasta_sequences = FASTAUtils.parse_file(
            expected_fasta_file_path
        )
        for fasta_sequence, expected_fasta_sequence in zip(
            fasta_sequences, expected_fasta_sequences
        ):
            assert fasta_sequence == expected_fasta_sequence

    @staticmethod
    def _validate_export_multi_sequence_fasta_file_from_fasta_sequences(
        sequence_fasta_file_path: pathlib.Path,
        export_directory_path: pathlib.Path,
    ):
        expected_fasta_sequences = FASTAUtils.parse_file(
            sequence_fasta_file_path
        )
        file_path = export_directory_path / 'exported-sequence.fasta'
        exported_fasta_file_path = FASTAUtils.export_as_fasta_file(
            file_path=file_path,
            sequences=expected_fasta_sequences,
        )
        TestFastaUtils._validate_multi_sequence_exported_fasta_file(
            exported_fasta_file_path=exported_fasta_file_path,
            expected_fasta_file_path=sequence_fasta_file_path,
        )

    @staticmethod
    def test_export_multi_sequence_fasta_file_from_fasta_dna_sequences(
        dna_sequences_fasta_file_path: pathlib.Path,
        tmp_path: pathlib.Path,
    ):
        TestFastaUtils._validate_export_multi_sequence_fasta_file_from_fasta_sequences(
            sequence_fasta_file_path=dna_sequences_fasta_file_path,
            export_directory_path=tmp_path,
        )

    @staticmethod
    def test_export_multi_sequence_fasta_file_from_fasta_rna_sequence(
        rna_sequences_fasta_file_path: pathlib.Path,
        tmp_path: pathlib.Path,
    ):
        TestFastaUtils._validate_export_multi_sequence_fasta_file_from_fasta_sequences(
            sequence_fasta_file_path=rna_sequences_fasta_file_path,
            export_directory_path=tmp_path,
        )

    @staticmethod
    def test_export_multi_sequence_fasta_file_from_fasta_protein_sequences(
        protein_sequences_fasta_file_path: pathlib.Path,
        tmp_path: pathlib.Path,
    ):
        TestFastaUtils._validate_export_multi_sequence_fasta_file_from_fasta_sequences(
            sequence_fasta_file_path=protein_sequences_fasta_file_path,
            export_directory_path=tmp_path,
        )
