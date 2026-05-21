import filecmp
import pathlib
import pytest

from pytest_datadir.plugin import LazyDataDir

from orfaqs.modules.python.orfaqsproteindiscovery.orfaqsproteindiscovery import (
    GenomicSequence,
    NucleotideUtils,
    ORFaqsProteinsDiscoveryApi,
    ORFaqsDiscoveredProteinRecord,
    RNAReadingFrame,
)


@pytest.fixture
def discovered_proteins_directory(
    lazy_shared_datadir: LazyDataDir,
) -> pathlib.Path:
    return lazy_shared_datadir.original_datadir / 'discovered-proteins'


@pytest.fixture
def discovered_proteins_samples_directory(
    discovered_proteins_directory: pathlib.Path,
) -> pathlib.Path:
    return discovered_proteins_directory / 'samples'


@pytest.fixture
def input_sequences_path(
    lazy_shared_datadir: LazyDataDir,
) -> pathlib.Path:
    return lazy_shared_datadir.original_datadir / 'input-sequences'


@pytest.fixture
def input_dna_sequence_str(
    input_sequences_path: pathlib.Path,
) -> str:
    file_path = input_sequences_path / 'dna-sequence.txt'
    with open(file_path, 'r', encoding='utf-8') as i_file:
        return i_file.read()

    return None


@pytest.fixture
def input_rna_sequence_str(
    input_sequences_path: pathlib.Path,
) -> str:
    file_path = input_sequences_path / 'rna-sequence.txt'
    with open(file_path, 'r', encoding='utf-8') as i_file:
        return i_file.read()

    return None


@pytest.fixture
def input_dna_sequence_fasta_file(
    input_sequences_path: pathlib.Path,
) -> pathlib.Path:
    return input_sequences_path / 'dna-sequence.fasta'


@pytest.fixture
def input_dna_sequence(
    input_dna_sequence_str: str,
) -> GenomicSequence:
    return NucleotideUtils.make_sequence_object(input_dna_sequence_str)


@pytest.fixture
def input_rna_sequence(
    input_rna_sequence_str: str,
) -> GenomicSequence:
    return NucleotideUtils.make_sequence_object(input_rna_sequence_str)


@pytest.fixture
def discovered_proteins_expected_results_dir(
    discovered_proteins_directory: pathlib.Path,
) -> pathlib.Path:
    return discovered_proteins_directory / 'expected-results'


class TestORFaqsProteinDiscoveryApi:
    @staticmethod
    def test_default_export_format():
        assert ORFaqsProteinsDiscoveryApi.default_export_format() == 'csv'

    @staticmethod
    def test_default_output_directory():
        assert ORFaqsProteinsDiscoveryApi.default_output_directory() == './'

    @staticmethod
    def test_available_export_formats():
        expected_export_formats = ['csv', 'json', 'xlsx']
        available_export_formats = sorted(
            ORFaqsProteinsDiscoveryApi.available_export_formats()
        )
        assert available_export_formats == expected_export_formats

    @staticmethod
    def test_exported_dataframe_keys():
        expected_exported_dataframe_keys: list[str] = [
            'index',
            'source_uid',
            'strand_type',
            'reading_frame',
            'genomic_sequence_position',
            'genomic_sequence',
            'protein',
            'protein_length',
        ]
        assert (
            expected_exported_dataframe_keys
            == ORFaqsProteinsDiscoveryApi.exported_dataframe_keys()
        )

    @staticmethod
    def test_find_discovered_protein_files(
        discovered_proteins_samples_directory: pathlib.Path,
    ):
        expected_discovered_protein_file_names: list[pathlib.Path] = [
            'discovered-proteins-0.csv',
            'discovered-proteins-1.csv',
            'discovered-proteins-2.csv',
        ]
        discovered_proteins_files = (
            ORFaqsProteinsDiscoveryApi.find_discovered_protein_files(
                discovered_proteins_samples_directory
            )
        )
        assert len(expected_discovered_protein_file_names) == len(
            discovered_proteins_files
        )
        assert all(
            discovered_proteins_file.name
            in expected_discovered_protein_file_names
            for discovered_proteins_file in discovered_proteins_files
        )

    @staticmethod
    def test_unreferenced_uid():
        assert (
            ORFaqsProteinsDiscoveryApi.unreferenced_uid()
            == 'unknown_reference'
        )

    @staticmethod
    def _validate_discovered_proteins(
        discovered_proteins: list[ORFaqsDiscoveredProteinRecord]
        | pathlib.Path,
        number_proteins: int,
        expected_results_file: pathlib.Path,
    ):
        if isinstance(discovered_proteins, pathlib.Path):
            filecmp.cmp(
                discovered_proteins,
                expected_results_file,
                shallow=False,
            )
        elif isinstance(discovered_proteins, list):
            expected_results: list[str] = []
            with open(expected_results_file, 'r', encoding='utf-8') as i_file:
                for line in i_file:
                    expected_results.append(line.strip('\r\n'))
            expected_number_proteins = len(expected_results)
            assert number_proteins == expected_number_proteins
            assert len(discovered_proteins) == expected_number_proteins
            for discovered_protein, expected_result_str in zip(
                discovered_proteins, expected_results
            ):
                # All genomic sequences must begin with ATG or AUG.
                first_triplet = discovered_protein.genomic_sequence.triplet(0)
                assert ('ATG' == first_triplet) or ('AUG' == first_triplet)
                # All proteins must begin with Methionine.
                assert 'M' == discovered_protein.protein[0]
                assert (
                    discovered_protein.condensed_record_json_str
                    == expected_result_str
                )

    @staticmethod
    def test_discover_proteins_from_dna_sequence(
        input_dna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_dna_sequence,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_dna_sequence_str(
        input_dna_sequence_str: str,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_dna_sequence_str,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_rna_sequence(
        input_rna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_rna_sequence,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_first_frame_only(
        input_rna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_rna_sequence,
                frames=[RNAReadingFrame.FIRST_FRAME],
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins-first-reading-frame-only.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_second_frame_only(
        input_rna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_rna_sequence,
                frames=[RNAReadingFrame.SECOND_FRAME],
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins-second-reading-frame-only.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_third_frame_only(
        input_rna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_rna_sequence,
                frames=[RNAReadingFrame.THIRD_FRAME],
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins-third-reading-frame-only.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_no_reverse_complement(
        input_rna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_rna_sequence,
                include_reverse_complement=False,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins-no-reverse-complement.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def _discover_proteins_from_sequence_export_results(
        input_sequence: GenomicSequence,
        output_directory: pathlib.Path,
        expected_results_file: pathlib.Path,
        export_format: str,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_sequence,
                export_results=True,
                output_directory=output_directory,
                export_format=export_format,
            )
        )
        assert discovered_proteins.parent == output_directory
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_export_default_results(
        input_rna_sequence: GenomicSequence,
        tmp_path: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.csv'
        )
        TestORFaqsProteinDiscoveryApi._discover_proteins_from_sequence_export_results(
            input_sequence=input_rna_sequence,
            output_directory=tmp_path,
            expected_results_file=expected_results_file,
            export_format=None,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_export_csv_results(
        input_rna_sequence: GenomicSequence,
        tmp_path: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        export_format = 'csv'
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / f'expected-discovered-proteins.{export_format}'
        )
        TestORFaqsProteinDiscoveryApi._discover_proteins_from_sequence_export_results(
            input_sequence=input_rna_sequence,
            output_directory=tmp_path,
            expected_results_file=expected_results_file,
            export_format=export_format,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_export_json_results(
        input_rna_sequence: GenomicSequence,
        tmp_path: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        export_format = 'json'
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / f'expected-discovered-proteins.{export_format}'
        )
        TestORFaqsProteinDiscoveryApi._discover_proteins_from_sequence_export_results(
            input_sequence=input_rna_sequence,
            output_directory=tmp_path,
            expected_results_file=expected_results_file,
            export_format=export_format,
        )

    @staticmethod
    def test_discover_proteins_from_sequence_export_xlsx_results(
        input_rna_sequence: GenomicSequence,
        tmp_path: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        export_format = 'xlsx'
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / f'expected-discovered-proteins.{export_format}'
        )
        TestORFaqsProteinDiscoveryApi._discover_proteins_from_sequence_export_results(
            input_sequence=input_rna_sequence,
            output_directory=tmp_path,
            expected_results_file=expected_results_file,
            export_format=export_format,
        )

    @staticmethod
    def test_discover_proteins_from_rna_sequence_use_cpu(
        input_rna_sequence: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_sequence(
                genomic_sequence=input_rna_sequence,
                enable_gpu=False,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_fasta_file(
        input_dna_sequence_fasta_file: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_fasta_file(
                genomic_sequence=input_dna_sequence_fasta_file
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_from_fasta_file_exported_results(
        input_dna_sequence_fasta_file: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins_from_fasta_file(
                genomic_sequence=input_dna_sequence_fasta_file,
                export_results=True,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.csv'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_detect_sequence_str(
        input_rna_sequence_str: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins(
                genomic_sequence=input_rna_sequence_str
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_detect_sequence_str_export_results(
        input_rna_sequence_str: GenomicSequence,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins(
                genomic_sequence=input_rna_sequence_str,
                export_results=True,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.csv'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_detect_fasta_file(
        input_dna_sequence_fasta_file: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins(
                genomic_sequence=input_dna_sequence_fasta_file
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.txt'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )

    @staticmethod
    def test_discover_proteins_detect_fasta_file_export_results(
        input_dna_sequence_fasta_file: pathlib.Path,
        discovered_proteins_expected_results_dir: pathlib.Path,
    ):
        (discovered_proteins, number_proteins) = (
            ORFaqsProteinsDiscoveryApi.discover_proteins(
                genomic_sequence=input_dna_sequence_fasta_file,
                export_results=True,
            )
        )
        expected_results_file = (
            discovered_proteins_expected_results_dir
            / 'expected-discovered-proteins.csv'
        )
        TestORFaqsProteinDiscoveryApi._validate_discovered_proteins(
            discovered_proteins,
            number_proteins,
            expected_results_file=expected_results_file,
        )
