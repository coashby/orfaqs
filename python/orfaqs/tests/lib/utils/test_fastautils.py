"""
Fasta Utils test module.
"""

from orfaqs.lib.utils.fastautils import (
    FASTASequenceType,
    FASTASequence,
    FASTAUtils,
)


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
    def test_read_rna_sequences(shared_datadir):
        fasta_file_path = shared_datadir / 'test_rna_sequences.fasta'
        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        TestFastaUtils._validate_number_sequences(
            fasta_file_path,
            fasta_sequences,
        )
        for fasta_sequence in fasta_sequences:
            assert fasta_sequence.sequence_type == FASTASequenceType.RNA

    @staticmethod
    def test_read_dna_sequences(shared_datadir):
        fasta_file_path = shared_datadir / 'test_dna_sequences.fasta'
        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        TestFastaUtils._validate_number_sequences(
            fasta_file_path,
            fasta_sequences,
        )
        for fasta_sequence in fasta_sequences:
            assert fasta_sequence.sequence_type == FASTASequenceType.DNA

    @staticmethod
    def test_read_protein_sequences(shared_datadir):
        fasta_file_path = shared_datadir / 'test_protein_sequences.fasta'
        fasta_sequences = FASTAUtils.parse_file(fasta_file_path)
        TestFastaUtils._validate_number_sequences(
            fasta_file_path,
            fasta_sequences,
        )
        for fasta_sequence in fasta_sequences:
            assert fasta_sequence.sequence_type == FASTASequenceType.AMINO_ACID

    @staticmethod
    def test_invalid_sequence(shared_datadir):
        fasta_file_path = shared_datadir / 'test_invalid_sequence.fasta'
        value_error_exception_thrown = False
        try:
            _ = FASTAUtils.parse_file(fasta_file_path)

        except ValueError:
            value_error_exception_thrown = True

        assert value_error_exception_thrown
