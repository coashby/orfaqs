"""
Ribosomes
"""

import enum
import logging
import numpy as np

import orfaqs.lib.core.aminoacids as _aminoacids
import orfaqs.lib.core.codons as _codons

from collections.abc import Iterator

from orfaqs.lib.utils.computeutils import (
    ComputeAccelerator,
    ComputeUtils,
)

from orfaqs.lib.core.aminoacids import AminoAcid
from orfaqs.lib.core.codons import Codon, CodonUtils
from orfaqs.lib.core.nucleotides import RNASequence
from orfaqs.lib.core.proteins import Protein

_logger = logging.getLogger(__name__)


class _ProfilerFunctionName(enum.Enum):
    RIBOSOME_TRANSLATE_RNA = enum.auto()
    RIBOSOME_UTILS_PROTEIN_DISCOVERY_ALL = enum.auto()
    RIBOSOME_UTILS_PROTEIN_DISCOVERY_FRAME = enum.auto()
    RIBOSOME_UTILS_FIND_START_CODONS = enum.auto()


_AMINO_ACID_ASSOCIATED_CODON_LUT: dict[AminoAcid, list[Codon]] = {
    _aminoacids.ALANINE: [_codons.GCU, _codons.GCC, _codons.GCA, _codons.GCG],
    _aminoacids.ARGININE: [
        _codons.CGU,
        _codons.CGC,
        _codons.CGA,
        _codons.CGG,
        _codons.AGA,
        _codons.AGG,
    ],
    _aminoacids.ASPARAGINE: [_codons.AAU, _codons.AAC],
    _aminoacids.ASPARTATE: [_codons.GAU, _codons.GAC],
    _aminoacids.CYSTEINE: [_codons.UGU, _codons.UGC],
    _aminoacids.GLUTAMINE: [_codons.CAA, _codons.CAG],
    _aminoacids.GLUTAMATE: [_codons.GAA, _codons.GAG],
    _aminoacids.GLYCINE: [_codons.GGU, _codons.GGC, _codons.GGA, _codons.GGG],
    _aminoacids.HISTIDINE: [_codons.CAU, _codons.CAC],
    _aminoacids.ISOLEUCINE: [_codons.AUU, _codons.AUC, _codons.AUA],
    _aminoacids.LEUCINE: [
        _codons.UUA,
        _codons.UUG,
        _codons.CUU,
        _codons.CUC,
        _codons.CUA,
        _codons.CUG,
    ],
    _aminoacids.LYSINE: [_codons.AAA, _codons.AAG],
    _aminoacids.METHIONINE: [_codons.AUG],
    _aminoacids.PHENYLALANINE: [_codons.UUU, _codons.UUC],
    _aminoacids.PROLINE: [_codons.CCU, _codons.CCC, _codons.CCA, _codons.CCG],
    _aminoacids.SERINE: [
        _codons.UCU,
        _codons.UCC,
        _codons.UCA,
        _codons.UCG,
        _codons.AGU,
        _codons.AGC,
    ],
    _aminoacids.THREONINE: [
        _codons.ACU,
        _codons.ACC,
        _codons.ACA,
        _codons.ACG,
    ],
    _aminoacids.TRYPTOPHAN: [_codons.UGG],
    _aminoacids.TYROSINE: [_codons.UAU, _codons.UAC],
    _aminoacids.VALINE: [_codons.GUU, _codons.GUC, _codons.GUA, _codons.GUG],
    _aminoacids.SELENOCYSTEINE: [_codons.UGA],
    _aminoacids.PYRROLYSINE: [_codons.UAG],
}


def _create_codon_to_amino_acid_lut() -> dict[Codon, AminoAcid]:
    lut: dict[Codon, AminoAcid] = {}
    for amino_acid, codons in _AMINO_ACID_ASSOCIATED_CODON_LUT.items():
        for codon in codons:
            lut[codon] = amino_acid

    return lut


_CODON_AMINO_ACID_LUT = _create_codon_to_amino_acid_lut()

_TRIPLET_AMINO_ACID_LUT = {
    codon.sequence_str: amino_acid
    for (codon, amino_acid) in _CODON_AMINO_ACID_LUT.items()
}


class RNAReadingFrame(enum.Enum):
    FIRST_FRAME = 1
    SECOND_FRAME = 2
    THIRD_FRAME = 3

    def __str__(self):
        return str(self.value)


_RNA_READING_FRAMES = [reading_frame for reading_frame in RNAReadingFrame]


class Ribosome:
    """Ribosome"""

    @staticmethod
    def translate_rna(
        rna_sequence: str | RNASequence,
        start_codons: list[str | Codon] = None,
        stop_codons: list[str | Codon] = None,
    ) -> Protein:
        """Translate an RNA sequence into an amino acid chain (protein)
        Note: In some cases, organisms or specific genome sites, such as
        mitochondrial genomes, may use different start and stop codons.

        Keyword Arguments:
        rna_sequence -- The RNA sequence to translate.
        start_codons -- A list of codons used by the ribosome to begin
                        translation. (Default: None)
        stop_codons -- A list of codons used by the ribosome to terminate
                       the RNA sequence translation process. (Default: None)
        """

        if not isinstance(start_codons, list):
            start_codons = RibosomeUtils.start_codon_triplets()
        if not isinstance(stop_codons, list):
            stop_codons = RibosomeUtils.stop_codon_triplets()

        protein: Protein = None
        codon_iterator = RibosomeUtils.read_codons(rna_sequence)
        for codon in codon_iterator:
            if codon in start_codons:
                protein = Protein(_CODON_AMINO_ACID_LUT.get(codon))
                break

        # Bind the remaining amino acids to the protein
        for codon in codon_iterator:
            if codon in stop_codons:
                # Expected termination of the sequence.
                protein.end_sequence()
                return protein

            protein.add_amino_acid(_CODON_AMINO_ACID_LUT.get(codon))

        message = (
            f'[INFO] {RibosomeUtils.__class__.__qualname__}: '
            'No stop codon was found during the formation of '
            'the protein.'
        )
        _logger.info(message)

        return protein


class RibosomeUtils:
    """RibosomeUtils"""

    @staticmethod
    def _compute_accelerator() -> ComputeAccelerator:
        compute_accelerator = ComputeUtils.get_default_compute_accelerator()
        if ComputeUtils.is_metal_compute_accelerator(compute_accelerator):
            compute_accelerator.load_kernel(
                'orfaqs/lib/core/kernels/metal/ribosomes.metal'
            )

        return compute_accelerator

    @staticmethod
    def start_codons() -> list[Codon]:
        return [_codons.AUG]

    @staticmethod
    def stop_codons() -> list[Codon]:
        return [_codons.UAA, _codons.UAG, _codons.UGA]

    @staticmethod
    def start_codon_triplets() -> list[str]:
        return [_codons.AUG.sequence_str]

    @staticmethod
    def stop_codon_triplets() -> list[str]:
        return [
            _codons.UAA.sequence_str,
            _codons.UAG.sequence_str,
            _codons.UGA.sequence_str,
        ]

    @staticmethod
    def available_reading_frames() -> list[RNAReadingFrame]:
        return _RNA_READING_FRAMES

    @staticmethod
    def codon_to_amino_acid(codon: Codon) -> AminoAcid:
        return _CODON_AMINO_ACID_LUT.get(codon)

    @staticmethod
    def triplet_to_amino_acid(triplet: str) -> AminoAcid:
        return _TRIPLET_AMINO_ACID_LUT.get(triplet)

    @staticmethod
    def sequence_start_stop_indices(
        rna_sequence: str | RNASequence, frame: RNAReadingFrame
    ) -> tuple[int, int]:
        sequence_length = len(rna_sequence)
        if sequence_length < CodonUtils.number_bases_per_codon():
            return None

        start_index = frame.value - 1
        stop_index = sequence_length - (
            (sequence_length - start_index)
            % CodonUtils.number_bases_per_codon()
        )

        return (start_index, stop_index)

    @staticmethod
    def _find_start_codons_gpu(
        rna_sequence: str | RNASequence,
        start_codons: list[Codon] = None,
    ) -> list[int]:
        kernel_function_name = 'find_start_codons'
        rna_sequence_str = rna_sequence
        if isinstance(rna_sequence, RNASequence):
            rna_sequence_str = rna_sequence.sequence_str

        start_condon_packed_str: str = ''
        for codon in start_codons:
            start_condon_packed_str += codon.sequence_str

        # Pad the codon list with dummy values
        # to ensure all codons are processed.
        base_padding_count = (
            3 - len(rna_sequence_str) % CodonUtils.number_bases_per_codon()
        ) % 3
        rna_sequence_str += 'x' * base_padding_count
        # start_codon_indices: list[bool] = [False] *
        number_codons = int(
            len(rna_sequence_str) / CodonUtils.number_bases_per_codon()
        )
        number_start_codons_buffer = [len(start_codons)]
        start_codon_indices: list[bool] = [False] * number_codons
        results_buffer_index = 3
        rna_sequence_uint8_list = list(rna_sequence_str.encode())

        #######################################################################
        # Create a new ComputeAccelerator object to process the RNA sequence.
        compute_accelerator = RibosomeUtils._compute_accelerator()
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=0,
            arg=rna_sequence_uint8_list,
            dtype=np.uint8,
        )
        start_condon_uint8_list = list(start_condon_packed_str.encode())
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=1,
            arg=start_condon_uint8_list,
            dtype=np.uint8,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=2,
            arg=number_start_codons_buffer,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            results_buffer_index,
            start_codon_indices,
            dtype=np.bool,
        )
        compute_accelerator.execute()
        start_codon_indices = compute_accelerator.get_buffer(
            buffer_index=results_buffer_index,
            dtype=np.bool,
        )

        start_codon_indices = np.where(start_codon_indices)[0].tolist()
        return start_codon_indices

    @staticmethod
    def _find_start_codons(
        rna_sequence: str | RNASequence,
        start_codons: list[Codon] = None,
    ) -> list[int]:
        start_codon_indices: list[int] = []
        codon_index = 0
        for codon in RibosomeUtils.read_codons(rna_sequence):
            if codon in start_codons:
                start_codon_indices.append(codon_index)

            codon_index += CodonUtils.number_bases_per_codon()

        return start_codon_indices

    @staticmethod
    def find_start_codons(
        rna_sequence: str | RNASequence,
        start_codons: list[Codon] = None,
        use_gpu: bool = True,
    ) -> list[int]:
        if start_codons is None:
            start_codons = RibosomeUtils.start_codons()

        if use_gpu:
            return RibosomeUtils._find_start_codons_gpu(
                rna_sequence,
                start_codons,
            )

        return RibosomeUtils._find_start_codons(
            rna_sequence,
            start_codons,
        )

    @staticmethod
    def read_triplets(rna_sequence: str | RNASequence) -> Iterator[str]:
        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)
        i = 0
        codon_step = CodonUtils.number_bases_per_codon()
        # Only read a multiple of codon_step bases.
        read_length = rna_sequence.sequence_length - (
            rna_sequence.sequence_length % codon_step
        )
        while i < read_length:
            yield str(rna_sequence[i : (i + codon_step)])
            i += codon_step

    @staticmethod
    def read_codons(rna_sequence: str | RNASequence) -> Iterator[Codon]:
        if rna_sequence is None or rna_sequence == '':
            return None

        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)
        i = 0
        codon_step = CodonUtils.number_bases_per_codon()
        # Only read a multiple of codon_step bases.
        read_length = rna_sequence.sequence_length - (
            rna_sequence.sequence_length % codon_step
        )
        while i < read_length:
            yield CodonUtils.base_triplet_to_codon(
                rna_sequence[i : (i + codon_step)]
            )
            i += codon_step
