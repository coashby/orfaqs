import orfaqs.lib.core.aminoacids as aminoacids
import pandas as pd


class TestAminoAcidsUtils:
    _amino_acid_utils = aminoacids.AminoAcidUtils()

    def test_available_amino_acids(self):
        reference_amino_acids_file_path = (
            './orfaqs/tests/resources/reference_amino_acids.csv'
        )
        reference_amino_acids_dataframe = pd.read_csv(
            reference_amino_acids_file_path
        )
        reference_amino_acids_dataframe = (
            reference_amino_acids_dataframe.set_index('name')
        )
        reference_amino_acids: dict[str, dict] = {}
        for amino_acid_name, row in reference_amino_acids_dataframe.iterrows():
            reference_amino_acids[amino_acid_name] = dict(row)

        available_amino_acids = self._amino_acid_utils.available_amino_acids()
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
