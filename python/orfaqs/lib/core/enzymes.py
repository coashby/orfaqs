"""
Enzymes
"""

from orfaqs.lib.core.nucleotides import (
    DNASequence,
    NucleicAcid,
    RNASequence,
    THYMINE,
    URACIL,
)


class RNAPolymerase:
    """RNAPolymerase"""

    @staticmethod
    def transcribe(
        dna_sequence: str | list[str] | list[NucleicAcid] | DNASequence,
        transcription_start_site: int = 0,
        termination_sequence: (
            str | list[str] | list[NucleicAcid] | DNASequence
        ) = None,
        rna_sequence_name: str = None,
    ) -> RNASequence:
        if not isinstance(dna_sequence, DNASequence):
            dna_sequence = DNASequence(dna_sequence)

        if (termination_sequence is not None) and (
            not isinstance(termination_sequence, DNASequence)
        ):
            termination_sequence = DNASequence(termination_sequence)

        dna_region: DNASequence = dna_sequence[transcription_start_site:]

        termination_site = dna_region.find_sequence(termination_sequence)
        if termination_site is None:
            termination_site = dna_region.sequence_length

        dna_region = dna_region[:termination_site]
        rna_sequence_str = dna_region.sequence_str.replace(
            THYMINE.symbol.lower(), URACIL.symbol.lower()
        )
        return RNASequence(rna_sequence_str, rna_sequence_name)
