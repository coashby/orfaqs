"""
Enzymes
"""

from orfaqs.libs.python.core.nucleotides import (
    DNASequence,
    StrandType,
    NucleicAcid,
    NucleotideUtils,
    RNASequence,
    THYMINE,
    URACIL,
)
from orfaqs.libs.python.core.sequence import SequenceUtils


class RNAPolymerase:
    """RNAPolymerase"""

    @staticmethod
    def transcribe(
        dna_sequence: str | list[str] | list[NucleicAcid] | DNASequence,
        strand_type: StrandType = None,
        transcription_start_site: int = 0,
        termination_sequence: (
            str | list[str] | list[NucleicAcid] | DNASequence
        ) = None,
        rna_sequence_name: str = None,
    ) -> RNASequence:
        if not isinstance(dna_sequence, DNASequence):
            dna_sequence = DNASequence(dna_sequence, strand_type=strand_type)

        if (termination_sequence is not None) and (
            not isinstance(termination_sequence, DNASequence)
        ):
            termination_sequence = DNASequence(
                termination_sequence, strand_type
            )

        dna_region: DNASequence = dna_sequence[transcription_start_site:]

        termination_site = NucleotideUtils.find_sequence(
            dna_region, termination_sequence
        )
        if termination_site is None:
            termination_site = dna_region.sequence_length

        dna_region = dna_region[:termination_site]
        rna_sequence_str = SequenceUtils.replace_symbols(
            dna_region,
            THYMINE.symbol,
            URACIL.symbol,
        )
        return RNASequence(
            rna_sequence_str,
            strand_type=strand_type,
            name=rna_sequence_name,
        )

    @staticmethod
    def reverse_transcribe(
        rna_sequence: str | list[str] | list[NucleicAcid] | RNASequence,
        strand_type: StrandType = None,
        transcription_start_site: int = 0,
        termination_sequence: (
            str | list[str] | list[NucleicAcid] | RNASequence
        ) = None,
        dna_sequence_name: str = None,
    ) -> RNASequence:
        if not isinstance(rna_sequence, RNASequence):
            rna_sequence = RNASequence(rna_sequence, strand_type=strand_type)

        if (termination_sequence is not None) and (
            not isinstance(termination_sequence, RNASequence)
        ):
            termination_sequence = RNASequence(
                termination_sequence, strand_type
            )

        rna_region: RNASequence = rna_sequence[transcription_start_site:]

        termination_site = NucleotideUtils.find_sequence(
            rna_region, termination_sequence
        )
        if termination_site is None:
            termination_site = rna_region.sequence_length

        rna_region = rna_region[:termination_site]
        dna_sequence_str = SequenceUtils.replace_symbols(
            rna_region,
            URACIL.symbol,
            THYMINE.symbol,
        )
        return DNASequence(
            dna_sequence_str,
            strand_type=strand_type,
            name=dna_sequence_name,
        )
