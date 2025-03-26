'''
Ribosomes
'''
import logging
import orfaqs.lib.core.aminoacids as _aminoacids
import orfaqs.lib.core.codons as _codons
from enum import Enum
from tqdm import tqdm

from orfaqs.lib.core.aminoacids import AminoAcid
from orfaqs.lib.core.codons import (
    Codon,
    CodonUtils
)
from orfaqs.lib.core.nucleotides import (
    NucleicAcid,
    RNASequence
)
from orfaqs.lib.core.proteins import Protein

_logger = logging.getLogger(__name__)

_AMINO_ACID_ASSOCIATED_CODON_LUT: dict[AminoAcid, list[Codon]] = {
    _aminoacids.ALANINE: [
        _codons.GCU,
        _codons.GCC,
        _codons.GCA,
        _codons.GCG
    ],
    _aminoacids.ARGININE: [
        _codons.CGU,
        _codons.CGC,
        _codons.CGA,
        _codons.CGG,
        _codons.AGA,
        _codons.AGG
    ],
    _aminoacids.ASPARAGINE: [
        _codons.AAU,
        _codons.AAC
    ],
    _aminoacids.ASPARTATE: [
        _codons.GAU,
        _codons.GAC
    ],
    _aminoacids.CYSTEINE: [
        _codons.UGU,
        _codons.UGC
    ],
    _aminoacids.GLUTAMINE: [
        _codons.CAA,
        _codons.CAG
    ],
    _aminoacids.GLUTAMATE: [
        _codons.GAA,
        _codons.GAG
    ],
    _aminoacids.GLYCINE: [
        _codons.GGU,
        _codons.GGC,
        _codons.GGA,
        _codons.GGG
    ],
    _aminoacids.HISTIDINE: [
        _codons.CAU,
        _codons.CAC
    ],
    _aminoacids.ISOLEUCINE: [
        _codons.AUU,
        _codons.AUC,
        _codons.AUA
    ],
    _aminoacids.LEUCINE: [
        _codons.UUA,
        _codons.UUG,
        _codons.CUU,
        _codons.CUC,
        _codons.CUA,
        _codons.CUG
    ],
    _aminoacids.LYSINE: [
        _codons.AAA,
        _codons.AAG
    ],
    _aminoacids.METHIONINE: [
        _codons.AUG
    ],
    _aminoacids.PHENYLALANINE: [
        _codons.UUU,
        _codons.UUC
    ],
    _aminoacids.PROLINE: [
        _codons.CCU,
        _codons.CCC,
        _codons.CCA,
        _codons.CCG
    ],
    _aminoacids.SERINE: [
        _codons.UCU,
        _codons.UCC,
        _codons.UCA,
        _codons.UCG,
        _codons.AGU,
        _codons.AGC
    ],
    _aminoacids.THREONINE: [
        _codons.ACU,
        _codons.ACC,
        _codons.ACA,
        _codons.ACG
    ],
    _aminoacids.TRYPTOPHAN: [
        _codons.UGG
    ],
    _aminoacids.TYROSINE: [
        _codons.UAU,
        _codons.UAC
    ],
    _aminoacids.VALINE: [
        _codons.GUU,
        _codons.GUC,
        _codons.GUA,
        _codons.GUG
    ],
    _aminoacids.SELENOCYSTEINE: [
        _codons.UGA
    ],
    _aminoacids.PYRROLYSINE: [
        _codons.UAG
    ]
}


def _create_codon_to_amino_acid_lut() -> dict[Codon, AminoAcid]:
    lut: dict[Codon, AminoAcid] = {}
    for (amino_acid, codons) in _AMINO_ACID_ASSOCIATED_CODON_LUT.items():
        for codon in codons:
            lut[codon] = amino_acid

    return lut


_CODON_AMINO_ACID_LUT = _create_codon_to_amino_acid_lut()


class RNAReadingFrame(Enum):
    FIRST_FRAME = 1
    SECOND_FRAME = 2
    THIRD_FRAME = 3


_RNA_READING_FRAMES = [reading_frame for reading_frame in RNAReadingFrame]


class Ribosome:
    '''Ribosome'''
    @staticmethod
    def translate_rna(rna_sequence: str | RNASequence,
                      start_codons: list[Codon] = None,
                      stop_codons: list[Codon] = None) -> Protein:
        '''Translate an RNA sequence into an amino acid chain (protein)
        Note: In some cases, organisms or specific genome sites, such as
        mitochondrial genomes, may use different start and stop codons.

        Keyword Arguments:
        rna_sequence -- The RNA sequence to translate.
        start_codons -- A list of codons used by the ribosome to begin
                        translation. (Default: None)
        stop_codons -- A list of codons used by the ribosome to terminate
                       the RNA sequence translation process. (Default: None)
        '''
        if not isinstance(start_codons, list):
            start_codons = RibosomeUtils.start_codons()
        if not isinstance(stop_codons, list):
            stop_codons = RibosomeUtils.stop_codons()

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

        message = (f'[INFO] {RibosomeUtils.__class__.__qualname__}: '
                   'No stop codon was found during the formation of '
                   'the protein.')
        _logger.info(message)
        return protein


class RibosomeUtils:
    '''RibosomeUtils'''

    @staticmethod
    def start_codons() -> list[Codon]:
        return [
            _codons.AUG
        ]

    @staticmethod
    def stop_codons() -> list[Codon]:
        return [
            _codons.UAA,
            _codons.UAG,
            _codons.UGA
        ]

    @staticmethod
    def available_reading_frames() -> list[RNAReadingFrame]:
        return _RNA_READING_FRAMES

    @staticmethod
    def codon_to_amino_acid(codon: str | Codon) -> AminoAcid:
        return _CODON_AMINO_ACID_LUT.get(codon)

    @staticmethod
    def sequence_start_stop_indices(
        rna_sequence: str | RNASequence,
            frame: RNAReadingFrame) -> tuple[int, int]:
        sequence_length = len(rna_sequence)
        if sequence_length < CodonUtils.number_bases_per_codon():
            return None

        start_index = frame.value - 1
        stop_index = (
            sequence_length -
            ((sequence_length - start_index) %
             CodonUtils.number_bases_per_codon())
        )

        return (start_index, stop_index)

    @staticmethod
    def find_start_codons(rna_sequence:  str | RNASequence,
                          start_codons: list[Codon] = None) -> list[int]:
        if start_codons is None:
            start_codons = RibosomeUtils.start_codons()

        start_codon_indices: list[int] = []
        codon_index = 0
        for codon in RibosomeUtils.read_codons(rna_sequence):
            if codon in start_codons:
                start_codon_indices.append(codon_index)

            codon_index += CodonUtils.number_bases_per_codon()

        return start_codon_indices

    @staticmethod
    def read_codons(rna_sequence: str | RNASequence):
        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)
        i = 0
        codon_step = CodonUtils.number_bases_per_codon()
        # Only read a multiple of codon_step bases.
        read_length = (rna_sequence.sequence_length -
                       (rna_sequence.sequence_length % codon_step))
        while i < read_length:
            yield CodonUtils.base_triplet_to_codon(
                rna_sequence[i:(i+codon_step)]
            )
            i += codon_step

    @staticmethod
    def discover_proteins(
            rna_sequence: str | list[str] | list[NucleicAcid] | RNASequence,
            frame: RNAReadingFrame = None,
            start_codons: list[Codon] = None,
            stop_codons: list[Codon] = None) -> tuple[dict, int]:
        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence)

        reading_frames: list[RNAReadingFrame] = []
        if frame is None:
            reading_frames = _RNA_READING_FRAMES
        else:
            reading_frames = [frame]

        protein_map = {reading_frame: [] for reading_frame in reading_frames}
        protein_count = 0
        for reading_frame in reading_frames:
            (start_index,
             stop_index) = RibosomeUtils.sequence_start_stop_indices(
                rna_sequence,
                reading_frame
            )
            rna_sequence_frame = rna_sequence[start_index:stop_index]
            # Get all start codon positions for the given reading frame.
            start_codon_indices = RibosomeUtils.find_start_codons(
                rna_sequence_frame,
                start_codons
            )

            # For every start codon found, translate the RNA region
            # beginning at that codon.
            process_description = ('Translating RNA from start codons in '
                                   f'reading frame {reading_frame.value}...')
            for start_codon_index in tqdm(start_codon_indices,
                                          desc=process_description):
                rna_coding_region = rna_sequence_frame[start_codon_index:]
                translated_protein = Ribosome.translate_rna(
                    rna_coding_region,
                    start_codons=start_codons,
                    stop_codons=stop_codons
                )
                if isinstance(translated_protein, Protein):
                    protein_map[reading_frame].append(translated_protein)
                    protein_count += 1

            print(f'Protein Count: {protein_count}')
        return (protein_map, protein_count)
