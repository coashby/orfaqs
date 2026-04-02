import pytest

from collections.abc import Iterator

import orfaqs.lib.core.codons as _codons

from orfaqs.lib.core.ribosomes import Ribosome, RibosomeUtils, RNAReadingFrame
from orfaqs.lib.core.codons import Codon
from orfaqs.lib.core.aminoacids import AminoAcid
from orfaqs.lib.core.nucleotides import RNASequence
from orfaqs.lib.core.proteins import Protein


@pytest.fixture
def default_rna_sequence_str():
    return 'AUGUGCAUGCUAA'


@pytest.fixture
def default_rna_sequence(default_rna_sequence_str: str):
    return RNASequence(default_rna_sequence_str)


class TestRibosome:
    """Test cases for Ribosome class"""

    @staticmethod
    def test_translate_rna_with_string_sequence(default_rna_sequence_str: str):
        """Test translation of RNA sequence provided as string"""

        protein = Ribosome.translate_rna(default_rna_sequence_str)
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_rna_sequence_object(
        default_rna_sequence: RNASequence,
    ):
        """Test translation of RNA sequence provided as RNASequence object"""

        protein = Ribosome.translate_rna(default_rna_sequence)
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_custom_start_codons(
        default_rna_sequence: RNASequence,
    ):
        """Test translation with custom start codons"""

        start_codons = ['AUG']
        protein = Ribosome.translate_rna(
            default_rna_sequence, start_codons=start_codons
        )
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_custom_stop_codons(
        default_rna_sequence: RNASequence,
    ):
        """Test translation with custom stop codons"""

        stop_codons = ['UAA']
        protein = Ribosome.translate_rna(
            default_rna_sequence, stop_codons=stop_codons
        )
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_both_custom_codons(
        default_rna_sequence: RNASequence,
    ):
        """Test translation with both custom start and stop codons"""

        start_codons = ['AUG']
        stop_codons = ['UAA']
        protein = Ribosome.translate_rna(
            default_rna_sequence,
            start_codons=start_codons,
            stop_codons=stop_codons,
        )
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_with_none_defaults(
        default_rna_sequence: RNASequence,
    ):
        """Test translation with default None parameters"""

        protein = Ribosome.translate_rna(default_rna_sequence)
        assert isinstance(protein, Protein)

    @staticmethod
    def test_translate_rna_empty_sequence():
        """Test translation of empty RNA sequence"""
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
        """Test start_codons returns correct codon list"""
        TestRibosomeUtils._validate_start_codons(
            RibosomeUtils.start_codons(),
            element_type=Codon,
        )

    @staticmethod
    def test_stop_codons():
        """Test stop_codons returns correct codon list"""
        TestRibosomeUtils._validate_stop_codons(
            RibosomeUtils.stop_codons(),
            element_type=Codon,
        )

    @staticmethod
    def test_start_codon_triplets():
        """Test start_codon_triplets returns string representations"""
        TestRibosomeUtils._validate_start_codons(
            RibosomeUtils.start_codon_triplets(),
            element_type=str,
        )

    @staticmethod
    def test_stop_codon_triplets():
        """Test stop_codon_triplets returns string representations"""
        TestRibosomeUtils._validate_stop_codons(
            RibosomeUtils.stop_codon_triplets(),
            element_type=str,
        )

    @staticmethod
    def test_available_reading_frames():
        """Test available_reading_frames returns all frames"""
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
        """Test codon_to_amino_acid with valid codon"""
        codon = _codons.AAA
        amino_acid = RibosomeUtils.codon_to_amino_acid(codon)
        assert isinstance(amino_acid, AminoAcid)

    @staticmethod
    def test_triplet_to_amino_acid_valid_triplet():
        """Test triplet_to_amino_acid with valid triplet string"""
        triplet = 'AUG'
        amino_acid = RibosomeUtils.triplet_to_amino_acid(triplet)
        assert isinstance(amino_acid, AminoAcid)

    @staticmethod
    def test_triplet_to_amino_acid_invalid_triplet():
        """Test triplet_to_amino_acid with invalid triplet"""
        triplet = 'XXX'
        amino_acid = RibosomeUtils.triplet_to_amino_acid(triplet)
        assert amino_acid is None

    @staticmethod
    def test_sequence_start_stop_indices_with_frame(
        default_rna_sequence: RNASequence,
    ):
        """Test sequence_start_stop_indices returns tuple of integers"""

        frame = RNAReadingFrame.FIRST_FRAME
        start_stop_indices = RibosomeUtils.sequence_start_stop_indices(
            default_rna_sequence, frame
        )
        assert isinstance(start_stop_indices, tuple)
        assert len(start_stop_indices) == 2
        assert all(isinstance(idx, int) for idx in start_stop_indices)

    @staticmethod
    def _create_rna_sequence_with_start_codons() -> tuple[
        RNASequence, list[int]
    ]:
        start_codon_indices = [
            index * Codon.number_bases() for index in [0, 3, 5]
        ]
        start_codon_indices = sorted(start_codon_indices)
        rna_sequence_length = (
            start_codon_indices[-1] + 2 * Codon.number_bases()
        )
        index = 0
        rna_sequence_str = ''
        while index < rna_sequence_length:
            if index in start_codon_indices:
                rna_sequence_str += RibosomeUtils.start_codons()[
                    0
                ].sequence_str
            else:
                rna_sequence_str += 'UUU'
            index += Codon.number_bases()
        return (RNASequence(rna_sequence_str), start_codon_indices)

    @staticmethod
    def _validate_found_start_codons(
        start_codon_indices: list[int],
        expected_start_codon_indices: list[int],
    ):
        assert isinstance(start_codon_indices, list)
        assert all(isinstance(idx, int) for idx in start_codon_indices)
        assert start_codon_indices == expected_start_codon_indices

    @staticmethod
    def test_find_start_codons_with_string_sequence():
        """Test find_start_codons with string RNA sequence"""

        (
            rna_sequence,
            expected_start_codon_indices,
        ) = TestRibosomeUtils._create_rna_sequence_with_start_codons()
        start_codon_indices = RibosomeUtils.find_start_codons(
            rna_sequence.sequence_str
        )
        TestRibosomeUtils._validate_found_start_codons(
            start_codon_indices,
            expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_with_rna_sequence_object():
        """Test find_start_codons with RNASequence object"""

        (
            rna_sequence,
            expected_start_codon_indices,
        ) = TestRibosomeUtils._create_rna_sequence_with_start_codons()
        start_codon_indices = RibosomeUtils.find_start_codons(rna_sequence)
        TestRibosomeUtils._validate_found_start_codons(
            start_codon_indices,
            expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_with_gpu():
        """Test find_start_codons with GPU acceleration enabled"""

        (
            rna_sequence,
            expected_start_codon_indices,
        ) = TestRibosomeUtils._create_rna_sequence_with_start_codons()
        start_codon_indices = RibosomeUtils.find_start_codons(
            rna_sequence,
            use_gpu=True,
        )
        TestRibosomeUtils._validate_found_start_codons(
            start_codon_indices,
            expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_without_gpu():
        """Test find_start_codons with GPU acceleration disabled"""

        (
            rna_sequence,
            expected_start_codon_indices,
        ) = TestRibosomeUtils._create_rna_sequence_with_start_codons()
        start_codon_indices = RibosomeUtils.find_start_codons(
            rna_sequence,
            use_gpu=False,
        )
        TestRibosomeUtils._validate_found_start_codons(
            start_codon_indices,
            expected_start_codon_indices,
        )

    @staticmethod
    def test_find_start_codons_custom_start_codons():
        """Test find_start_codons with custom start codons"""

        custom_codons = [_codons.AUG]

        (
            rna_sequence,
            expected_start_codon_indices,
        ) = TestRibosomeUtils._create_rna_sequence_with_start_codons()
        start_codon_indices = RibosomeUtils.find_start_codons(
            rna_sequence,
            start_codons=custom_codons,
        )
        TestRibosomeUtils._validate_found_start_codons(
            start_codon_indices,
            expected_start_codon_indices,
        )

    @staticmethod
    def _validate_read_triplets(triplet_iterator: Iterator[str]):
        assert hasattr(triplet_iterator, '__iter__')
        triplets = list(triplet_iterator)
        assert all(
            (isinstance(triplet, str) and len(triplet) == Codon.number_bases())
            for triplet in triplets
        )

    @staticmethod
    def test_read_triplets_with_string_sequence(
        default_rna_sequence_str: RNASequence,
    ):
        """Test read_triplets returns iterator of strings"""

        triplet_iterator = RibosomeUtils.read_triplets(
            default_rna_sequence_str
        )
        TestRibosomeUtils._validate_read_triplets(triplet_iterator)

    @staticmethod
    def test_read_triplets_with_rna_sequence_object(
        default_rna_sequence: RNASequence,
    ):
        """Test read_triplets with RNASequence object"""

        triplet_iterator = RibosomeUtils.read_triplets(default_rna_sequence)
        TestRibosomeUtils._validate_read_triplets(triplet_iterator)

    @staticmethod
    def _validate_read_codons(codon_iterator: Iterator[Codon]):
        assert hasattr(codon_iterator, '__iter__')
        triplets = list(codon_iterator)
        assert all(isinstance(triplet, Codon) for triplet in triplets)

    @staticmethod
    def test_read_codons_with_string_sequence(
        default_rna_sequence_str: str,
    ):
        """Test read_codons returns iterator of Codons"""

        codons_iterator = RibosomeUtils.read_codons(default_rna_sequence_str)
        TestRibosomeUtils._validate_read_codons(codons_iterator)

    @staticmethod
    def test_read_codons_with_rna_sequence_object(
        default_rna_sequence: RNASequence,
    ):
        """Test read_codons with RNASequence object"""

        codons_iterator = RibosomeUtils.read_codons(default_rna_sequence)
        TestRibosomeUtils._validate_read_codons(codons_iterator)
