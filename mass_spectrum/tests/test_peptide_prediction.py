""" Peptide Prediction module tests."""


from .constants import TEST_FASTA_FILE, TEST_MZML_FILE
from ms_package.peptide_prediction import PeptideSearch
import pytest
import pandas as pd

test = PeptideSearch(fasta_path=str(TEST_FASTA_FILE), mzml_path=str(TEST_MZML_FILE))


class TestPeptideSearch:
    """A test class which conducts unit tests for PeptideSearch class."""

    def test_peptide_search(self):
        """Checks whether protein_ids and peptide_ids are returned after SimpleSearchEngineAlgorithm."""

        assert test.fasta_path.endswith('.fasta') is True
        assert test.mzml_path.endswith('.mzML') is True
        protein_ids = test.peptide_search()[0]
        peptide_ids = test.peptide_search()[1]
        assert isinstance(protein_ids, list)
        assert isinstance(peptide_ids, list)

    def test_get_peptide_identification_values(self):
        """Checks whether the PeptideSearch class has correct peptide values stored in dictionary."""

        peptide_ids = test.peptide_search()[1]
        dict_peptide_info = test.get_peptide_identification_values(peptide_ids=peptide_ids)
        assert isinstance(dict_peptide_info, dict)
        #  tests for first dictionary i.e. first peptide seq identified
        assert dict_peptide_info[0]['Peptide ID m/z'] == float(443.71)
        assert dict_peptide_info[0]['Peptide ID rt'] == float(1738.03)
        assert dict_peptide_info[0]['Peptide hit sequence'] == str('DDSPDLPK')
        assert dict_peptide_info[0]['Peptide hit score'] == float(0.03)

    def test_get_sequence(self):
        """Checks whether list of hit peptide sequences is correctly identified by algorithm"""

        peptide_ids = test.peptide_search()[1]
        sequences = test.get_sequence(peptide_ids=peptide_ids)
        assert isinstance(sequences, list)
        #  example of peptide sequences from test data only
        assert len(sequences) == 15
        assert sequences == ['DDSPDLPK', 'YIC(Carbamidomethyl)DNQDTISSK', 'C(Carbamidomethyl)C(Carbamidomethyl)TESLVNR',
                             'LC(Carbamidomethyl)VLHEK', 'DLGEEHFK', 'LC(Carbamidomethyl)VLHEK', 'LVTDLTK',
                             'DLGEEHFK', 'AEFVEVTK', 'EAC(Carbamidomethyl)FAVEGPK', 'EAC(Carbamidomethyl)FAVEGPK',
                             'GAC(Carbamidomethyl)LLPK', 'YLYEIAR', 'LVVSTQTALA', 'YLYEIAR']

    def test_peptide_wrapper(self):
        """Checks the peptide values stored in dataframe."""

        info_dataframe = test.peptide_wrapper()[0]
        assert isinstance(info_dataframe, pd.DataFrame)

        # examples of random values retrieved from dataframe
        info_dataframe = test.peptide_wrapper()[0]
        row_0 = next(info_dataframe.iterrows())[1]  # checks for values stored in 1st row i.e values for first peptide seq
        assert row_0['Peptide ID m/z'] == float(443.71)
        assert row_0['Peptide ID rt'] == float(1738.03)
        assert row_0['Peptide hit sequence'] == str('DDSPDLPK')
        assert row_0['Peptide hit score'] == float(0.03)
