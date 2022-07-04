"""Reader module tests."""

import pytest
import pandas as pd
from ms_package.reader import Reader
import xml.dom.minidom
import argparse
from .constants import TEST_FASTA_FILE, TEST_MZML_FILE, TEST_MZXML_FILE


test1 = Reader(str(TEST_MZML_FILE))
test2 = Reader(str(TEST_MZXML_FILE))
test3 = Reader(str(TEST_FASTA_FILE))


class TestReader:
    """A test class which conducts pytests on the Reader class."""

    def test_check_extension(self):
        """Tests whether the check_extension function returns True if a valid file format is passed and
        False if an invalid format is passed. Checks if the self.format values are set correctly."""
        assert test1.check_extension() is True
        assert test1.format == 'mzml'
        assert test2.check_extension() is True
        assert test2.format == 'mzxml'
        assert test3.check_extension() is False

    def test_parse_file(self):
        """Tests whether the parse_file method parses the input file and returns an xml.dom.minidom.Document object."""
        file1 = test1.parse_file()
        assert isinstance(file1, xml.dom.minidom.Document)
        file2 = test2.parse_file()
        assert isinstance(file2, xml.dom.minidom.Document)
        with pytest.raises(argparse.ArgumentTypeError):
            test3.parse_file()

    def test_get_spectrum_list(self):
        """Tests whether the get_spectrum_list method correctly extracts the spectra and returns an
        xml.dom.minidom.Element object."""
        file1 = test1.parse_file()
        output1 = test1.get_spectrum_list(file1)
        assert isinstance(output1, list)
        assert isinstance(output1[0], xml.dom.minidom.Element)
        assert len(output1) == 1684
        file2 = test2.parse_file()
        output2 = test2.get_spectrum_list(file2)
        assert isinstance(output2, list)
        assert isinstance(output2[0], xml.dom.minidom.Element)
        assert len(output2) == 7161

    def test_get_spectrum_dict(self):
        """Tests whether the get_spectrum_dict method correctly creates a dictionary of the spectra and returns an
        xml.dom.minidom.Element object."""
        file1 = test1.parse_file()
        list1 = test1.get_spectrum_list(file1)
        dict1 = test1.get_spectrum_dict(list1)
        assert isinstance(dict1, dict)
        assert len(dict1) == 1684
        assert isinstance(list(dict1.keys())[0], int)
        assert isinstance(dict1[0], xml.dom.minidom.Element)
        file2 = test2.parse_file()
        list2 = test2.get_spectrum_list(file2)
        dict2 = test2.get_spectrum_dict(list2)
        assert isinstance(dict2, dict)
        assert len(dict2) == 7161
        assert isinstance(list(dict2.keys())[0], int)
        assert isinstance(dict2[1], xml.dom.minidom.Element)

    def test_get_compression(self):
        """Tests whether the get_compression method extracts the compression from the input file correctly and creates
        the compression dictionary correctly. Checks whether an error is raised if run on an mzXML file."""
        file1 = test1.parse_file()
        list1 = test1.get_spectrum_list(file1)
        dict1 = test1.get_spectrum_dict(list1)
        test1.get_compression(dict1)
        assert isinstance(test1.compression, dict)
        assert list(test1.compression[0].keys()) == ['mz', 'intensity']
        assert len(test1.compression) == 1684
        assert test1.compression[0]['mz'] == {'data_type': '64-bit float', 'compression': 'no compression'}
        file2 = test2.parse_file()
        list2 = test2.get_spectrum_list(file2)
        dict2 = test2.get_spectrum_dict(list2)
        with pytest.raises(argparse.ArgumentTypeError):
            test2.get_compression(dict2)

    def test_get_binary_spectrum_values(self):
        """Tests whether the get_binary_spectrum_values method extracts the binary arrays from the input file and
        sets the self.binary_values attribute correctly. Checks whether an error is raised if run on an mzXML file."""
        file1 = test1.parse_file()
        list1 = test1.get_spectrum_list(file1)
        dict1 = test1.get_spectrum_dict(list1)
        test1.get_binary_spectrum_values(dict1)
        assert isinstance(test1.binary_values, dict)
        assert list(test1.binary_values[0].keys()) == ['mz', 'intensity']
        assert len(test1.binary_values) == 1684
        file2 = test2.parse_file()
        list2 = test2.get_spectrum_list(file2)
        dict2 = test2.get_spectrum_dict(list2)
        with pytest.raises(argparse.ArgumentTypeError):
            test2.get_binary_spectrum_values(dict2)

    def test_get_values(self):
        """Tests whether the base peak m/z, base peak intensity, total ion current and the lowest and highest observed
        m/z values are extracted correctly."""
        file1 = test1.parse_file()
        list1 = test1.get_spectrum_list(file1)
        dict1 = test1.get_spectrum_dict(list1)
        values1 = test1.get_values(dict1)
        assert isinstance(values1, dict)
        assert len(values1[0]) == 6
        keys = ['spectra_id', 'base_peak_m/z', 'base_peak_intensity', 'total_ion_current', 'lowest_observed_m/z',
                'highest_observed_m/z']
        assert keys == list(values1[0].keys())
        id0 = {'spectra_id': 0,
               'base_peak_m/z': 391.28,
               'base_peak_intensity': 928844.25,
               'total_ion_current': 6937649.0,
               'lowest_observed_m/z': 300.0,
               'highest_observed_m/z': 2008.46}
        assert values1[0] == id0
        assert test1.values == values1

    def test_analyse_spectrum(self):
        """Tests whether the wrapper method analyse_spectrum returns a pandas dataframe."""
        result1 = test1.analyse_spectrum()
        assert isinstance(result1, pd.DataFrame)
