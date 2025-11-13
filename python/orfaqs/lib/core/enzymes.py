"""
Enzymes
"""

from orfaqs.lib.core.nucleotides import (
    DNASequence,
    StrandType,
    NucleicAcid,
    NucleotideUtils,
    RNASequence,
    THYMINE,
    URACIL,
)


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
        rna_sequence_str = dna_region.sequence_str.replace(
            THYMINE.symbol.lower(), URACIL.symbol.lower()
        )
        return RNASequence(
            rna_sequence_str, strand_type=strand_type, name=rna_sequence_name
        )
