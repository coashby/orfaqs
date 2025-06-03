import orfaqs.lib.core.aminoacids as aminoacids
import pandas as pd
import pytest


@pytest.fixture
def reference_amino_acids(shared_datadir) -> dict[str, dict]:
    reference_amino_acids_file_path = (
        shared_datadir / 'reference_amino_acids.csv'
    )
    reference_amino_acids_dataframe = pd.read_csv(
        reference_amino_acids_file_path
    )

    reference_amino_acids_dataframe = (
        reference_amino_acids_dataframe.set_index('name')
    )
    amino_acids: dict[str, dict] = {}
    for amino_acid_name, row in reference_amino_acids_dataframe.iterrows():
        amino_acids[amino_acid_name] = dict(row)

    return amino_acids


class TestAminoAcidsUtils:
    @staticmethod
    def test_available_amino_acids(reference_amino_acids: dict[str, dict]):
        available_amino_acids = (
            aminoacids.AminoAcidUtils.available_amino_acids()
        )
        assert len(reference_amino_acids) == len(available_amino_acids)
        for amino_acid in available_amino_acids:
            reference_abbreviations = reference_amino_acids.get(
                amino_acid.name
            )
            # Assert all properties are the same.
            assert reference_abbreviations is not None
            assert (
                amino_acid.abbreviation
                == reference_abbreviations['abbreviation']
            )
            assert amino_acid.symbol == reference_abbreviations['symbol']
