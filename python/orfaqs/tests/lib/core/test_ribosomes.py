import json
import pathlib
import pytest

from collections.abc import Iterator
from pytest_datadir.plugin import LazyDataDir

import orfaqs.lib.core.codons as _codons

from orfaqs.lib.core.ribosomes import Ribosome, RibosomeUtils, RNAReadingFrame
from orfaqs.lib.core.codons import Codon
from orfaqs.lib.core.aminoacids import AminoAcid
from orfaqs.lib.core.nucleotides import RNASequence
from orfaqs.lib.core.proteins import Protein


def _read_rna_sequence_txt_file(
    rna_sequence_file_path: str | pathlib.Path,
) -> str:
    rna_sequence_str: str = None
    with open(rna_sequence_file_path, 'r', encoding='utf-8') as i_file:
        rna_sequence_str = i_file.read()
        rna_sequence_str = rna_sequence_str.replace('\n', '')
        rna_sequence_str = rna_sequence_str.replace('\r', '')
        rna_sequence_str = rna_sequence_str.replace(' ', '')

    return rna_sequence_str


def _read_json_file(file_path: str | pathlib.Path) -> any:
    json_contents: any = None
    with open(file_path, 'r', encoding='utf-8') as i_file:
        json_contents = json.load(i_file)

    return json_contents


@pytest.fixture
def sample_rna_sequence_short_str(shared_datadir: LazyDataDir) -> str:
    rna_sequence_file_path = (
        shared_datadir
        / 'sample-rna-sequences'
        / 'sample-rna-sequence-short.txt'
    )
    return _read_rna_sequence_txt_file(rna_sequence_file_path)


@pytest.fixture
def sample_rna_sequence_short(
    sample_rna_sequence_short_str: str,
) -> RNASequence:
    return RNASequence(sample_rna_sequence_short_str)


@pytest.fixture
def sample_rna_sequence_long(shared_datadir: LazyDataDir) -> RNASequence:
    rna_sequence_file_path = (
        shared_datadir
        / 'sample-rna-sequences'
        / 'sample-rna-sequence-long.txt'
    )
    return RNASequence(_read_rna_sequence_txt_file(rna_sequence_file_path))


@pytest.fixture
def sample_rna_sequence_expected_results_path(
    shared_datadir: LazyDataDir,
) -> pathlib.Path:
    return shared_datadir / 'expected-results' / 'sample-rna-sequences'


@pytest.fixture
def sample_rna_sequence_short_expected_start_codon_indices(
    sample_rna_sequence_expected_results_path: pathlib.Path,
) -> list[int]:
    file_path = (
        sample_rna_sequence_expected_results_path
        / 'start-codon-indices-short.json'
    )
    return _read_json_file(file_path)


@pytest.fixture
def sample_rna_sequence_long_expected_start_codon_indices(
    sample_rna_sequence_expected_results_path: pathlib.Path,
) -> list[int]:
    file_path = (
        sample_rna_sequence_expected_results_path
        / 'start-codon-indices-long.json'
    )
    return _read_json_file(file_path)


@pytest.fixture
def sample_rna_sequence_short_expected_stop_codon_indices(
    sample_rna_sequence_expected_results_path: pathlib.Path,
) -> list[int]:
    file_path = (
        sample_rna_sequence_expected_results_path
        / 'stop-codon-indices-short.json'
    )
    return _read_json_file(file_path)


@pytest.fixture
def sample_rna_sequence_long_expected_stop_codon_indices(
    sample_rna_sequence_expected_results_path: pathlib.Path,
) -> list[int]:
    file_path = (
        sample_rna_sequence_expected_results_path
        / 'stop-codon-indices-long.json'
    )
    return _read_json_file(file_path)


def _format_expected_protein_results(results_object: dict):
    expected_orf_proteins: dict[int, Protein] = {}
    for index, protein_str in results_object.items():
        expected_orf_proteins[int(index)] = Protein(protein_str)

    return expected_orf_proteins


@pytest.fixture
def sample_rna_sequence_short_expected_orf_proteins(
    sample_rna_sequence_expected_results_path: pathlib.Path,
) -> dict[int, Protein]:
    file_path = (
        sample_rna_sequence_expected_results_path
        / 'all-orf-proteins-short.json'
    )
    results_object: dict = _read_json_file(file_path)
    return _format_expected_protein_results(results_object)


@pytest.fixture
def sample_rna_sequence_long_expected_orf_proteins(
    sample_rna_sequence_expected_results_path: pathlib.Path,
) -> dict[int, Protein]:
    file_path = (
        sample_rna_sequence_expected_results_path
        / 'all-orf-proteins-long.json'
    )
    results_object: dict = _read_json_file(file_path)
    return _format_expected_protein_results(results_object)


class TestRibosome:
    """Test cases for Ribosome class"""

    @staticmethod
    def test_translate_rna_with_sequence_str_input(
        sample_rna_sequence_short_str: str,
    ):
        """Test translation of RNA sequence provided as string."""

        protein = Ribosome.translate_rna(sample_rna_sequence_short_str)
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_sequence_object_input(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test translation of RNA sequence provided as RNASequence object."""

        protein = Ribosome.translate_rna(sample_rna_sequence_short)
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_custom_start_codons(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test translation with custom start codons."""

        start_codons = ['AUG']
        protein = Ribosome.translate_rna(
            sample_rna_sequence_short, start_codons=start_codons
        )
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_custom_stop_codons(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test translation with custom stop codons."""

        stop_codons = ['UAA']
        protein = Ribosome.translate_rna(
            sample_rna_sequence_short, stop_codons=stop_codons
        )
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_both_custom_codons(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test translation with both custom start and stop codons."""

        start_codons = ['AUG']
        stop_codons = ['UAA']
        protein = Ribosome.translate_rna(
            sample_rna_sequence_short,
            start_codons=start_codons,
            stop_codons=stop_codons,
        )
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_none_defaults(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test translation with default None parameters."""

        protein = Ribosome.translate_rna(sample_rna_sequence_short)
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_empty_sequence():
        """Test translation of empty RNA sequence."""
        rna_sequence_str = ''
        protein = Ribosome.translate_rna(rna_sequence_str)
        assert protein is None


class TestRibosomeUtils:
    """Test cases for RibosomeUtils class"""

    @staticmethod
    def _validate_start_codons(
        start_codons: list[Codon | str],
        element_type: type,
    ):
        assert isinstance(start_codons, list)
        assert 'AUG' in start_codons
        assert all(isinstance(codon, element_type) for codon in start_codons)

    @staticmethod
    def _validate_stop_codons(
        stop_codons: list[Codon | str],
        element_type: type,
    ):
        assert isinstance(stop_codons, list)
        for stop_codon in ['UAA', 'UAG', 'UGA']:
            assert stop_codon in stop_codons
        assert all(isinstance(codon, element_type) for codon in stop_codons)

    @staticmethod
    def test_start_codons():
        """Test start_codons returns correct codon list."""
        TestRibosomeUtils._validate_start_codons(
            RibosomeUtils.start_codons(),
            element_type=Codon,
        )

    @staticmethod
    def test_stop_codons():
        """Test stop_codons returns correct codon list."""
        TestRibosomeUtils._validate_stop_codons(
            RibosomeUtils.stop_codons(),
            element_type=Codon,
        )

    @staticmethod
    def test_start_codon_triplets():
        """Test start_codon_triplets returns string representations."""
        TestRibosomeUtils._validate_start_codons(
            RibosomeUtils.start_codon_triplets(),
            element_type=str,
        )

    @staticmethod
    def test_stop_codon_triplets():
        """Test stop_codon_triplets returns string representations."""
        TestRibosomeUtils._validate_stop_codons(
            RibosomeUtils.stop_codon_triplets(),
            element_type=str,
        )

    @staticmethod
    def test_available_reading_frames():
        """Test available_reading_frames returns all frames."""
        reading_frames = RibosomeUtils.available_reading_frames()
        assert isinstance(reading_frames, list)
        assert len(reading_frames) == 3
        assert all(
            isinstance(frame, RNAReadingFrame) for frame in reading_frames
        )
        assert RNAReadingFrame.FIRST_FRAME in reading_frames
        assert RNAReadingFrame.SECOND_FRAME in reading_frames
        assert RNAReadingFrame.THIRD_FRAME in reading_frames

    @staticmethod
    def test_codon_to_amino_acid_valid_codon():
        """Test codon_to_amino_acid with valid codon."""
        codon = _codons.AAA
        amino_acid = RibosomeUtils.codon_to_amino_acid(codon)
        assert isinstance(amino_acid, AminoAcid)

    @staticmethod
    def test_triplet_to_amino_acid_valid_triplet():
        """Test triplet_to_amino_acid with valid triplet string."""
        triplet = 'AUG'
        amino_acid = RibosomeUtils.triplet_to_amino_acid(triplet)
        assert isinstance(amino_acid, AminoAcid)

    @staticmethod
    def test_triplet_to_amino_acid_invalid_triplet():
        """Test triplet_to_amino_acid with invalid triplet."""
        triplet = 'XXX'
        amino_acid = RibosomeUtils.triplet_to_amino_acid(triplet)
        assert amino_acid is None

    @staticmethod
    def test_sequence_start_stop_indices_with_frame(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test sequence_start_stop_indices returns tuple of integers."""

        frame = RNAReadingFrame.FIRST_FRAME
        start_stop_indices = RibosomeUtils.sequence_start_stop_indices(
            sample_rna_sequence_short, frame
        )
        assert isinstance(start_stop_indices, tuple)
        assert len(start_stop_indices) == 2
        assert all(isinstance(idx, int) for idx in start_stop_indices)

    @staticmethod
    def _validate_read_triplets(triplet_iterator: Iterator[str]):
        assert hasattr(triplet_iterator, '__iter__')
        triplets = list(triplet_iterator)
        assert all(
            (isinstance(triplet, str) and len(triplet) == Codon.number_bases())
            for triplet in triplets
        )

    @staticmethod
    def test_read_triplets_using_squence_str_input(
        sample_rna_sequence_short_str: RNASequence,
    ):
        """Test read_triplets returns iterator of strings."""

        triplet_iterator = RibosomeUtils.read_triplets(
            sample_rna_sequence_short_str
        )
        TestRibosomeUtils._validate_read_triplets(triplet_iterator)

    @staticmethod
    def test_read_triplets_with_sequence_object_input(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test read_triplets with RNASequence object."""

        triplet_iterator = RibosomeUtils.read_triplets(
            sample_rna_sequence_short
        )
        TestRibosomeUtils._validate_read_triplets(triplet_iterator)

    @staticmethod
    def _validate_read_codons(codon_iterator: Iterator[Codon]):
        assert hasattr(codon_iterator, '__iter__')
        triplets = list(codon_iterator)
        assert all(isinstance(triplet, Codon) for triplet in triplets)

    @staticmethod
    def test_read_codons_with_sequence_str_input(
        sample_rna_sequence_short_str: str,
    ):
        """Test read_codons returns iterator of Codons."""

        codons_iterator = RibosomeUtils.read_codons(
            sample_rna_sequence_short_str
        )
        TestRibosomeUtils._validate_read_codons(codons_iterator)

    @staticmethod
    def test_read_codons_with_sequence_object_input(
        sample_rna_sequence_short: RNASequence,
    ):
        """Test read_codons with RNASequence object."""

        codons_iterator = RibosomeUtils.read_codons(sample_rna_sequence_short)
        TestRibosomeUtils._validate_read_codons(codons_iterator)

    @staticmethod
    def _validate_found_codons(
        codon_indices: list[int],
        expected_codon_indices: list[int],
    ):
        assert isinstance(codon_indices, list)
        assert all(isinstance(idx, int) for idx in codon_indices)
        assert codon_indices == expected_codon_indices

    @staticmethod
    def test_find_start_codons_with_sequence_str_input(
        sample_rna_sequence_short_str: str,
        sample_rna_sequence_short_expected_start_codon_indices: list[int],
    ):
        """Test find_start_codons with string RNA sequence."""

        start_codon_indices = RibosomeUtils.find_start_codons(
            sample_rna_sequence_short_str
        )
        TestRibosomeUtils._validate_found_codons(
            start_codon_indices,
            sample_rna_sequence_short_expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_with_sequence_object_input(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_start_codon_indices: list[int],
    ):
        """Test find_start_codons with RNASequence object."""

        start_codon_indices = RibosomeUtils.find_start_codons(
            sample_rna_sequence_short
        )
        TestRibosomeUtils._validate_found_codons(
            start_codon_indices,
            sample_rna_sequence_short_expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_with_gpu_enabled(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_start_codon_indices: list[int],
    ):
        """Test find_start_codons with GPU acceleration enabled."""

        start_codon_indices = RibosomeUtils.find_start_codons(
            sample_rna_sequence_short,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_found_codons(
            start_codon_indices,
            sample_rna_sequence_short_expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_long_sequence(
        sample_rna_sequence_long: RNASequence,
        sample_rna_sequence_long_expected_start_codon_indices: list[int],
    ):
        """Test find_start_codons with GPU acceleration disabled."""

        start_codon_indices = RibosomeUtils.find_start_codons(
            sample_rna_sequence_long,
            use_gpu=False,
        )
        TestRibosomeUtils._validate_found_codons(
            start_codon_indices,
            sample_rna_sequence_long_expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_with_gpu_enabled_long_sequence(
        sample_rna_sequence_long: RNASequence,
        sample_rna_sequence_long_expected_start_codon_indices: list[int],
    ):
        """Test find_start_codons with GPU acceleration disabled."""

        start_codon_indices = RibosomeUtils.find_start_codons(
            sample_rna_sequence_long,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_found_codons(
            start_codon_indices,
            sample_rna_sequence_long_expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_custom_start_codons(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_start_codon_indices: list[int],
    ):
        """Test find_start_codons with custom start codons."""

        custom_codons = [_codons.AUG]
        start_codon_indices = RibosomeUtils.find_start_codons(
            sample_rna_sequence_short,
            start_codons=custom_codons,
        )
        TestRibosomeUtils._validate_found_codons(
            start_codon_indices,
            sample_rna_sequence_short_expected_start_codon_indices,
        )

    @staticmethod
    def test_find_stop_codons_with_sequence_str_input(
        sample_rna_sequence_short_str: str,
        sample_rna_sequence_short_expected_stop_codon_indices: list[int],
    ):
        """Test find_stop_codons with string RNA sequence."""

        stop_codon_indices = RibosomeUtils.find_stop_codons(
            sample_rna_sequence_short_str
        )
        TestRibosomeUtils._validate_found_codons(
            stop_codon_indices,
            sample_rna_sequence_short_expected_stop_codon_indices,
        )

    @staticmethod
    def test_find_stop_codons_with_sequence_object_input(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_stop_codon_indices: list[int],
    ):
        """Test find_stop_codons with RNASequence object."""

        stop_codon_indices = RibosomeUtils.find_stop_codons(
            sample_rna_sequence_short
        )
        TestRibosomeUtils._validate_found_codons(
            stop_codon_indices,
            sample_rna_sequence_short_expected_stop_codon_indices,
        )

    @staticmethod
    def test_find_stop_codons_with_gpu_enabled(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_stop_codon_indices: list[int],
    ):
        """Test find_stop_codons with GPU acceleration enabled."""

        stop_codon_indices = RibosomeUtils.find_stop_codons(
            sample_rna_sequence_short,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_found_codons(
            stop_codon_indices,
            sample_rna_sequence_short_expected_stop_codon_indices,
        )

    @staticmethod
    def test_find_stop_codons_long_sequence(
        sample_rna_sequence_long: RNASequence,
        sample_rna_sequence_long_expected_stop_codon_indices: list[int],
    ):
        """Test find_stop_codons with GPU acceleration disabled."""

        stop_codon_indices = RibosomeUtils.find_stop_codons(
            sample_rna_sequence_long,
            use_gpu=False,
        )
        TestRibosomeUtils._validate_found_codons(
            stop_codon_indices,
            sample_rna_sequence_long_expected_stop_codon_indices,
        )

    @staticmethod
    def test_find_stop_codons_with_gpu_enabled_long_sequence(
        sample_rna_sequence_long: RNASequence,
        sample_rna_sequence_long_expected_stop_codon_indices: list[int],
    ):
        """Test find_stop_codons with GPU acceleration disabled."""

        stop_codon_indices = RibosomeUtils.find_stop_codons(
            sample_rna_sequence_long,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_found_codons(
            stop_codon_indices,
            sample_rna_sequence_long_expected_stop_codon_indices,
        )

    @staticmethod
    def test_find_stop_codons_custom_stop_codons(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_stop_codon_indices: list[int],
    ):
        """Test find_stop_codons with custom stop codons."""

        custom_codons = [_codons.UAA, _codons.UAG, _codons.UGA]
        stop_codon_indices = RibosomeUtils.find_stop_codons(
            sample_rna_sequence_short,
            stop_codons=custom_codons,
        )
        TestRibosomeUtils._validate_found_codons(
            stop_codon_indices,
            sample_rna_sequence_short_expected_stop_codon_indices,
        )

    @staticmethod
    def _validate_translate_all_orfs(
        all_orf_proteins: dict[int, Protein],
        expected_orf_proteins: dict[int, Protein],
    ):
        assert len(all_orf_proteins) == len(expected_orf_proteins)
        miss_count = 0
        print()
        miss_dict = {}
        for index, protein in all_orf_proteins.items():
            assert isinstance(protein, Protein)
            assert index in expected_orf_proteins
            if protein != expected_orf_proteins[index]:
                miss_count += 1
                miss_dict[index] = [
                    protein.sequence_str,
                    expected_orf_proteins[index].sequence_str,
                ]

        if miss_count > 0:
            print()

    @staticmethod
    def test_translate_all_orfs_with_sequence_str_input(
        sample_rna_sequence_short_str,
        sample_rna_sequence_short_expected_orf_proteins: dict[int, Protein],
    ):
        """Test retranslate_all_orfs with string RNASequence."""

        all_orf_proteins = RibosomeUtils.translate_all_orfs(
            rna_sequence=sample_rna_sequence_short_str,
            use_gpu=False,
        )
        TestRibosomeUtils._validate_translate_all_orfs(
            all_orf_proteins,
            sample_rna_sequence_short_expected_orf_proteins,
        )

    @staticmethod
    def test_translate_all_orfs_with_sequence_object_input(
        sample_rna_sequence_short,
        sample_rna_sequence_short_expected_orf_proteins: dict[int, Protein],
    ):
        """Test retranslate_all_orfs with RNASequence object."""

        all_orf_proteins = RibosomeUtils.translate_all_orfs(
            rna_sequence=sample_rna_sequence_short,
            use_gpu=False,
        )
        TestRibosomeUtils._validate_translate_all_orfs(
            all_orf_proteins,
            sample_rna_sequence_short_expected_orf_proteins,
        )

    @staticmethod
    def test_translate_all_orfs_with_gpu_enabled(
        sample_rna_sequence_short: RNASequence,
        sample_rna_sequence_short_expected_orf_proteins: dict[int, Protein],
    ):
        """Test retranslate_all_orfs with RNASequence object."""
        all_orf_proteins = RibosomeUtils.translate_all_orfs(
            rna_sequence=sample_rna_sequence_short,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_translate_all_orfs(
            all_orf_proteins,
            sample_rna_sequence_short_expected_orf_proteins,
        )

    @staticmethod
    def test_translate_all_orfs_long_sequence(
        sample_rna_sequence_long: RNASequence,
        sample_rna_sequence_long_expected_orf_proteins: dict[int, Protein],
    ):
        """Test retranslate_all_orfs with RNASequence object."""
        all_orf_proteins = RibosomeUtils.translate_all_orfs(
            rna_sequence=sample_rna_sequence_long,
            use_gpu=False,
        )
        TestRibosomeUtils._validate_translate_all_orfs(
            all_orf_proteins,
            sample_rna_sequence_long_expected_orf_proteins,
        )

    @staticmethod
    def test_translate_all_orfs_with_gpu_enabled_long_sequence(
        sample_rna_sequence_long: RNASequence,
        sample_rna_sequence_long_expected_orf_proteins: dict[int, Protein],
    ):
        """Test retranslate_all_orfs with RNASequence object."""
        all_orf_proteins = RibosomeUtils.translate_all_orfs(
            rna_sequence=sample_rna_sequence_long,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_translate_all_orfs(
            all_orf_proteins,
            sample_rna_sequence_long_expected_orf_proteins,
        )
