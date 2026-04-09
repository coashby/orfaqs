"""
Ribosomes
"""

import enum
import logging
import numpy as np

import orfaqs.lib.core.aminoacids as _aminoacids
import orfaqs.lib.core.codons as _codons

from collections.abc import Iterator
from tqdm import tqdm

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
        if sequence_length < Codon.number_bases():
            return None

        start_index = frame.value - 1
        stop_index = sequence_length - (
            (sequence_length - start_index) % Codon.number_bases()
        )

        return (start_index, stop_index)

    @staticmethod
    def _prepare_rna_sequence_for_gpu(
        rna_sequence: str | RNASequence,
    ) -> list[int]:
        rna_sequence_str = rna_sequence
        if isinstance(rna_sequence, RNASequence):
            rna_sequence_str = rna_sequence.sequence_str

        # Truncate the RNA sequence so that only valid codons are processed.
        sequence_length = len(rna_sequence_str)
        processing_length = sequence_length - (
            sequence_length % Codon.number_bases()
        )
        rna_sequence_str = rna_sequence_str[:processing_length]

        return list(rna_sequence_str.encode())

    @staticmethod
    def _prepare_codons_for_gpu(codons: list[Codon]) -> str:
        codons_packed_str: str = ''
        for codon in codons:
            codons_packed_str += codon.sequence_str

        return list(codons_packed_str.encode())

    @staticmethod
    def _find_codons_gpu(
        rna_sequence: str | RNASequence,
        codons: list[Codon],
    ) -> list[int]:
        rna_sequence_buffer = RibosomeUtils._prepare_rna_sequence_for_gpu(
            rna_sequence
        )
        reference_codons_buffer = RibosomeUtils._prepare_codons_for_gpu(codons)

        number_codons = int(len(rna_sequence_buffer) / Codon.number_bases())
        number_reference_codons = len(codons)
        found_codon_indices: list[bool] = [False] * number_codons
        results_buffer_index = 5

        #######################################################################
        # Create a new ComputeAccelerator object to process the RNA sequence.
        arg_index = ComputeAccelerator.arg_index_generator()
        compute_accelerator = RibosomeUtils._compute_accelerator()
        data_stride = compute_accelerator.calculate_data_stride(number_codons)
        kernel_function_name = 'find_codons'
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=rna_sequence_buffer,
            dtype=np.uint8,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=number_codons,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=data_stride,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=reference_codons_buffer,
            dtype=np.uint8,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=number_reference_codons,
            dtype=np.uint32,
        )
        results_buffer_index = next(arg_index)
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=results_buffer_index,
            arg=found_codon_indices,
            dtype=np.bool,
        )
        compute_accelerator.execute()
        found_codon_indices = compute_accelerator.get_buffer(
            buffer_index=results_buffer_index,
            dtype=np.bool,
        )

        found_codon_indices = [
            index * Codon.number_bases()
            for index in np.where(found_codon_indices)[0].tolist()
        ]
        return found_codon_indices

    @staticmethod
    def _find_codons(
        rna_sequence: str | RNASequence,
        codons: list[Codon],
    ) -> list[int]:
        found_codon_indices: list[int] = []
        codon_index = 0
        number_codons = int(len(rna_sequence) / Codon.number_bases())
        description = 'Finding codons...'
        with tqdm(total=number_codons, desc=description) as progress_bar:
            for codon in RibosomeUtils.read_codons(rna_sequence):
                if codon in codons:
                    found_codon_indices.append(codon_index)

                codon_index += Codon.number_bases()
                progress_bar.update(1)

        return found_codon_indices

    @staticmethod
    def find_start_codons(
        rna_sequence: str | RNASequence,
        start_codons: list[Codon] = None,
        use_gpu: bool = False,
    ) -> list[int]:
        """
        Returns a sorted list of base indices from the input sequence, in
        ascending order, indicating the position of start codons found within
        the sequence.
        ----------
        Arguments:
        ----------
        rna_sequence (str | RNASequence):
            The input RNA sequence to inspect for start codons.

        start_codons (list[Codon]):
            (Optional) A list of codons known to be, or labeled as, start
            codons.

        use_gpu (bool):
            (Optional) If True, then the GPU is used for finding the start
            codon indices. Otherwise, the CPU is used. The default value is
            False.
        """
        if start_codons is None:
            start_codons = RibosomeUtils.start_codons()

        if ComputeUtils.compute_accelerator_available() and use_gpu:
            return RibosomeUtils._find_codons_gpu(
                rna_sequence,
                start_codons,
            )

        return RibosomeUtils._find_codons(
            rna_sequence,
            start_codons,
        )

    @staticmethod
    def find_stop_codons(
        rna_sequence: str | RNASequence,
        stop_codons: list[Codon] = None,
        use_gpu: bool = False,
    ) -> list[int]:
        """
        Returns a sorted list of base indices from the input sequence, in
        ascending order, indicating the position of stop codons found within
        the sequence.
        ----------
        Arguments:
        ----------
        rna_sequence (str | RNASequence):
            The input RNA sequence to inspect for stop codons.

        stop_codons (list[Codon]):
            (Optional) A list of codons known to be, or labeled as, stop
            codons.

        use_gpu (bool):
            (Optional) If True, then the GPU is used for finding the stop
            codon indices. Otherwise, the CPU is used. The default value is
            False.
        """
        if stop_codons is None:
            stop_codons = RibosomeUtils.stop_codons()

        if ComputeUtils.compute_accelerator_available() and use_gpu:
            return RibosomeUtils._find_codons_gpu(
                rna_sequence,
                stop_codons,
            )

        return RibosomeUtils._find_codons(
            rna_sequence,
            stop_codons,
        )

    @staticmethod
    def read_triplets(rna_sequence: str | RNASequence) -> Iterator[str]:
        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)
        i = 0
        codon_step = Codon.number_bases()
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
        codon_step = Codon.number_bases()
        # Only read a multiple of codon_step bases.
        read_length = rna_sequence.sequence_length - (
            rna_sequence.sequence_length % codon_step
        )
        while i < read_length:
            yield CodonUtils.base_triplet_to_codon(
                rna_sequence[i : (i + codon_step)]
            )
            i += codon_step

    @staticmethod
    def _all_orf_lengths_gpu(
        start_codon_indices: list[int],
        stop_codon_indices: list[int],
    ) -> list[int]:
        number_start_codons = len(start_codon_indices)
        number_stop_codons = len(stop_codon_indices)
        orf_lengths: list[int] = [0] * number_start_codons
        #######################################################################
        # Create a new ComputeAccelerator object to process the indices.
        arg_index = ComputeAccelerator.arg_index_generator()
        compute_accelerator = RibosomeUtils._compute_accelerator()
        data_stride = compute_accelerator.calculate_data_stride(
            number_start_codons
        )
        kernel_function_name = 'all_orf_lengths'
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=start_codon_indices,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=number_start_codons,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=stop_codon_indices,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=number_stop_codons,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=data_stride,
            dtype=np.uint32,
        )
        results_buffer_index = next(arg_index)
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=results_buffer_index,
            arg=orf_lengths,
            dtype=np.int32,
        )

        compute_accelerator.execute()
        orf_lengths = compute_accelerator.get_buffer(
            buffer_index=results_buffer_index,
            dtype=np.int32,
        )
        return orf_lengths

    @staticmethod
    def _create_uint8_buffer(buffer_length: int) -> list[int]:
        return list((' ' * buffer_length).encode())

    @staticmethod
    def _rna_to_amino_acid_sequence_buffer_gpu(
        rna_sequence: str | RNASequence,
    ) -> list[int]:
        rna_sequence_buffer = RibosomeUtils._prepare_rna_sequence_for_gpu(
            rna_sequence
        )
        number_codons = int(len(rna_sequence_buffer) / Codon.number_bases())
        amino_acid_sequence_length = number_codons
        amino_acid_sequence_buffer = RibosomeUtils._create_uint8_buffer(
            amino_acid_sequence_length
        )

        #######################################################################
        # Create a new ComputeAccelerator object to process the
        # amino acid sequence.
        arg_index = ComputeAccelerator.arg_index_generator()
        compute_accelerator = RibosomeUtils._compute_accelerator()
        # Use the largest data dimension available to calculate the
        # data stride.
        data_stride = compute_accelerator.calculate_data_stride(number_codons)
        kernel_function_name = 'rna_to_amino_acid_sequence'
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=rna_sequence_buffer,
            dtype=np.uint8,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=number_codons,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=data_stride,
            dtype=np.uint32,
        )
        results_buffer_index = next(arg_index)
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=results_buffer_index,
            arg=amino_acid_sequence_buffer,
            dtype=np.uint8,
        )

        compute_accelerator.execute()
        amino_acid_sequence_buffer = compute_accelerator.get_buffer(
            buffer_index=results_buffer_index,
            dtype=np.uint8,
        )

        return amino_acid_sequence_buffer

    @staticmethod
    def _translate_all_orfs_gpu(
        rna_sequence: str | RNASequence,
        start_codon_indices: list[int],
        stop_codon_indices: list[int],
    ) -> dict[int, Protein]:
        amino_acid_sequence_buffer = (
            RibosomeUtils._rna_to_amino_acid_sequence_buffer_gpu(rna_sequence)
        )
        # Sort and convert all codon indices to condon-centric indices.
        start_codon_indices = RibosomeUtils._convert_base_to_codon_indices(
            sorted(start_codon_indices)
        )
        stop_codon_indices = RibosomeUtils._convert_base_to_codon_indices(
            sorted(stop_codon_indices)
        )

        # Compute the lengths of each ORF.
        orf_lengths = RibosomeUtils._all_orf_lengths_gpu(
            start_codon_indices,
            stop_codon_indices,
        )
        # Remove all start codon indices that do not start a valid ORF. This is
        # determined by finding the first start codon index that resulted in an
        # ORF length of 0.
        if 0 in orf_lengths:
            end_index = orf_lengths.index(0)
            orf_lengths = orf_lengths[:end_index]
            start_codon_indices = start_codon_indices[:end_index]

        number_start_codons = len(start_codon_indices)
        # Create the protein buffer.
        protein_buffer_length = sum(orf_lengths)
        all_orf_proteins_buffer = RibosomeUtils._create_uint8_buffer(
            protein_buffer_length
        )
        #######################################################################
        # Create a new ComputeAccelerator object to process the
        # amino acid sequence.
        arg_index = ComputeAccelerator.arg_index_generator()
        compute_accelerator = RibosomeUtils._compute_accelerator()
        # Use the largest data dimension available to calculate the
        # data stride.
        data_stride = compute_accelerator.calculate_data_stride(
            number_start_codons
        )
        kernel_function_name = 'translate_all_orfs'
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=amino_acid_sequence_buffer,
            dtype=np.uint8,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=start_codon_indices,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=orf_lengths,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=number_start_codons,
            dtype=np.uint32,
        )
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=next(arg_index),
            arg=data_stride,
            dtype=np.uint32,
        )
        results_buffer_index = next(arg_index)
        compute_accelerator.set_arg(
            kernel_function_name,
            arg_index=results_buffer_index,
            arg=all_orf_proteins_buffer,
            dtype=np.uint8,
        )

        compute_accelerator.execute()
        all_orf_proteins_buffer = compute_accelerator.get_buffer(
            buffer_index=results_buffer_index,
            dtype=np.uint8,
        )
        all_orf_proteins: dict[int, Protein] = {}
        buffer_start_index = 0

        for index, length in enumerate(orf_lengths):
            buffer_end_index = buffer_start_index + length
            protein_buffer = all_orf_proteins_buffer[
                buffer_start_index:buffer_end_index
            ]
            protein_str = bytearray(protein_buffer).decode().replace(' ', '')
            base_index = start_codon_indices[index] * Codon.number_bases()
            all_orf_proteins[base_index] = Protein(protein_str)
            buffer_start_index = buffer_end_index
        return all_orf_proteins

    @staticmethod
    def _convert_base_to_codon_indices(codon_indices: list[int]) -> list[int]:
        return [int(index / Codon.number_bases()) for index in codon_indices]

    @staticmethod
    def _translate_all_orfs(
        rna_sequence: str | RNASequence,
        start_codon_indices: list[int],
        stop_codon_indices: list[int],
    ) -> dict[int, Protein]:
        # Sort and convert all codon indices to condon-centric indices.
        start_codon_indices = RibosomeUtils._convert_base_to_codon_indices(
            sorted(start_codon_indices)
        )
        stop_codon_indices = RibosomeUtils._convert_base_to_codon_indices(
            sorted(stop_codon_indices)
        )

        # Convert all codons to their corresponding amino acid
        amino_acids: list[AminoAcid] = [
            RibosomeUtils.codon_to_amino_acid(codon)
            for codon in RibosomeUtils.read_codons(rna_sequence)
        ]

        # Copy all valid amino acid sequences into the ORF protein list.
        all_orf_proteins: dict[int, Protein] = {}
        # Mark the end of the translation routine by adding a stop index to
        # the start_codon_indices list.
        end_of_list_marker = -1
        start_codon_indices.append(end_of_list_marker)
        # Reverse the indices to provide stack functionality.
        start_codon_indices.reverse()
        stop_codon_indices.reverse()
        start_codon_index = start_codon_indices.pop()
        stop_codon_index = stop_codon_indices.pop()
        number_start_codons = len(start_codon_indices)
        description = 'Translating ORFs...'
        with tqdm(total=number_start_codons, desc=description) as progress_bar:
            while start_codon_index != end_of_list_marker:
                if start_codon_index > stop_codon_index:
                    if len(stop_codon_indices) == 0:
                        # There are no more stop codons to mark the
                        # end of a valid reading frame.
                        break
                    else:
                        # Get a new stop codon index.
                        stop_codon_index = stop_codon_indices.pop()
                        # Reevaluate outer loop conditions before
                        # attempting transcription.
                        continue

                # Transcribe the ORF.
                start_codon_base_index = (
                    start_codon_index * Codon.number_bases()
                )
                # Do not include the stop codon in the translation.
                all_orf_proteins[start_codon_base_index] = Protein(
                    amino_acids[start_codon_index:stop_codon_index]
                )
                # Get a new start codon index.
                start_codon_index = start_codon_indices.pop()
                progress_bar.update(1)

        return all_orf_proteins

    @staticmethod
    def translate_all_orfs(
        rna_sequence: str | RNASequence,
        start_codons: list[Codon] = None,
        stop_codons: list[Codon] = None,
        use_gpu: bool = False,
    ) -> dict[int, Protein]:
        """
        Returns a dictionary of RNA sequence base indices mapped to the protein
        encoded by the corresponding ORF.
        ----------
        Arguments:
        ----------
        rna_sequence (str | RNASequence):
            The input RNA sequence to inspect for stop codons.

        start_codons (list[Codon]):
            (Optional) A list of codons known to be, or labeled as, start
            codons.

        stop_codons (list[Codon]):
            (Optional) A list of codons known to be, or labeled as, stop
            codons.

        use_gpu (bool):
            (Optional) If True, then the GPU is used for finding the stop
            codon indices. Otherwise, the CPU is used. The default value is
            False.
        """
        if start_codons is None:
            start_codons = RibosomeUtils.start_codons()
        if stop_codons is None:
            stop_codons = RibosomeUtils.stop_codons()

        start_codon_indices = RibosomeUtils.find_start_codons(
            rna_sequence=rna_sequence,
            start_codons=start_codons,
            use_gpu=use_gpu,
        )
        stop_codon_indices = RibosomeUtils.find_stop_codons(
            rna_sequence=rna_sequence,
            stop_codons=stop_codons,
            use_gpu=use_gpu,
        )

        if ComputeUtils.compute_accelerator_available() and use_gpu:
            return RibosomeUtils._translate_all_orfs_gpu(
                rna_sequence=rna_sequence,
                start_codon_indices=start_codon_indices,
                stop_codon_indices=stop_codon_indices,
            )

        return RibosomeUtils._translate_all_orfs(
            rna_sequence=rna_sequence,
            start_codon_indices=start_codon_indices,
            stop_codon_indices=stop_codon_indices,
        )
