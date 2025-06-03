import orfaqs.lib.core.codons as codons
import pandas as pd
import pytest


@pytest.fixture
def reference_codons(shared_datadir) -> dict[str, dict]:
    reference_codons_file_path = shared_datadir / 'reference_codons.csv'
    reference_codons_dataframe = pd.read_csv(reference_codons_file_path)

    reference_codons_dataframe = reference_codons_dataframe.set_index(
        'triplet_code'
    )
    codons: dict[str, dict] = {}
    for triplet, row in reference_codons_dataframe.iterrows():
        codons[triplet] = dict(row)

    return codons


class TestCodonUtils:
    @staticmethod
    def test_available_codons(reference_codons: dict[str, dict]):
        available_codons: list[codons.Codon] = (
            codons.CodonUtils.available_codons()
        )
        assert len(available_codons) == len(reference_codons)
        for codon in available_codons:
            assert codon.sequence_str in reference_codons

    def test_base_triplet_to_codon(self, reference_codons: dict[str, dict]):
        for triplet_str in reference_codons.keys():
            codon = codons.CodonUtils.base_triplet_to_codon(triplet_str)
            assert isinstance(codon, codons.Codon)
