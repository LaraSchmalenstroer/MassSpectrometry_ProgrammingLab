""" Protein Prediction module tests. """

import os
import json

import pandas as pd
import pytest

from ms_package.startup import DATA_DIR
from ms_package.protein_prediction import ProteinSearch
from ms_package.peptide_prediction import PeptideSearch
from .constants import TEST_FASTA_FILE, TEST_MZML_FILE

# Test peptide list from test file
pep = PeptideSearch(fasta_path=str(TEST_FASTA_FILE), mzml_path=str(TEST_MZML_FILE))
peptide_list = pep.peptide_wrapper()[1]


class TestProteinSearch:
    """Unit tests for ProteinSearch class """

    def test_filter_peptides(self):
        """Checks whether peptide list is filtered
        """
        pro = ProteinSearch(peptide_list=peptide_list)
        filtered_pep_list = pro.filter_peptides()
        assert isinstance(filtered_pep_list, list)
        assert len(filtered_pep_list) == 8

    def test_divide_into_chunks(self):
        """Checks whether a given list is divided into list of 15 elements
        """
        pro = ProteinSearch(peptide_list=peptide_list)
        # example list only for this method
        list_for_chunks = ['DDSPDLPK', 'DLGEEHFK', 'LVTDLTK', 'DLGEEHFK', 'LVTDLTK', 'LVVSTQTALA', 'YLYEIAR',
                           'AEFVEVTK', 'YLYEIAR', 'LVVSTQTALA', 'YLYEIAR', 'DLGEEHFK', 'AEFVEVTK', 'YLYEIAR',
                           'DDSPDLPK', 'DLGEEHFK', 'YLYEIAR', 'LVVSTQTALA', 'YLYEIAR']
        pro.divide_into_chunks(list_for_chunks)

        assert len(pro.sel_peptides) == 2
        assert pro.sel_peptides == [['DDSPDLPK', 'DLGEEHFK', 'LVTDLTK', 'DLGEEHFK', 'LVTDLTK', 'LVVSTQTALA', 'YLYEIAR',
                                     'AEFVEVTK', 'YLYEIAR', 'LVVSTQTALA', 'YLYEIAR', 'DLGEEHFK', 'AEFVEVTK', 'YLYEIAR',
                                     'DDSPDLPK'], ['DLGEEHFK', 'YLYEIAR', 'LVVSTQTALA', 'YLYEIAR']]

    def test_get_api_query(self):
        """Checks if the API query is in the correct format
        """
        pro = ProteinSearch(peptide_list=peptide_list)
        filtered_pep_list = pro.filter_peptides()
        query = pro.get_api_query(filtered_pep_list[0:3])
        ans = 'DDSPDLPK%2CDLGEEHFK%2CLVTDLTK'
        assert query == ans

    def test_file_handle(self):
        """Checks if a data file is created"""
        pro = ProteinSearch(peptide_list=peptide_list)
        pro.file_handle()
        assert os.path.exists(os.path.join(DATA_DIR, "proteins", "protein_data.csv")) is True

    def test_proteins_api(self):
        """ Checks the response of The Proteins API
        """
        pro = ProteinSearch(peptide_list=peptide_list)
        pro.divide_into_chunks(["DLGEEHFK"])
        assert pro.sel_peptides == [["DLGEEHFK"]]

        pro.proteins_api()
        pro_response = pro.protein_response
        assert isinstance(pro_response, str)

        pro_response_json = json.loads(pro_response)[0]

        assert pro_response_json["features"][0]["peptide"] == "DLGEEHFK"
        assert pro_response_json["accession"] == "A0A140T897"
        assert pro_response_json["taxid"] == 9913
        assert f"{pro_response_json['features'][0]['begin']}-{pro_response_json['features'][0]['end']}" == "37-44"

    def test_get_proteins(self):
        """Check the wrapper function
        """
        # test full exists condition
        pro = ProteinSearch(peptide_list=["DLGEEHFK"])
        pro.get_proteins()
        assert isinstance(pro.ans_df, pd.DataFrame)
        assert pro.ans_df["Peptide"].unique() == ["DLGEEHFK"]

        # test full API condition
        os.remove(os.path.join(DATA_DIR, "proteins", "protein_data.csv"))  # delete cache database
        pro = ProteinSearch(peptide_list=["DLGEEHFK"])
        pro.get_proteins()
        assert isinstance(pro.ans_df, pd.DataFrame)
        assert pro.ans_df["Peptide"].unique() == ["DLGEEHFK"]

        # test half exists condition
        pro = ProteinSearch(peptide_list=["DLGEEHFK", "AEFVEVTK"])
        pro.get_proteins()
        assert isinstance(pro.ans_df, pd.DataFrame)
        assert list(pro.ans_df["Peptide"].unique()) == ["DLGEEHFK", "AEFVEVTK"]
