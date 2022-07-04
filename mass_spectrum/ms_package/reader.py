import xml.dom.minidom
from xml.dom import minidom as md
from typing import List, Dict
import base64
import zlib
import struct
import logging
import pandas as pd
import argparse


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Reader:
    """Parses the input mzml/mzXml file and extracts spectrum values."""
    def __init__(self, path):
        self.path = path  # path to input file
        self.format = None  # can be 'mzml' or 'mzxml', set in check_extension
        self.compression = None  # contains compression dict for each spectrum
        self.binary_values = None  # contains binary values (m/z and intensity arrays) for each spectrum id
        self.spectrum_data = None  # decoded intensity and m/z array values
        self.values = None  # base peak m/z, base peak intensity, lowest and highest observed m/z and total ion current

    def check_extension(self) -> bool:
        """Checks if the extension of the parsed file is either .mzML or .mzXML.

        Returns
        -------
        bool: True if extension is allowed

        """
        if self.path.endswith('.mzML'):
            self.format = 'mzml'
            return True
        elif self.path.endswith('.mzXML'):
            self.format = 'mzxml'
            return True
        else:
            return False

    def parse_file(self) -> xml.dom.minidom.Document:
        """Parses the input file and creates minidom object.

        Returns
        -------
        parsed_file: xml.dom.minidom.Document
            xml.dom.minidom.Document from the parsed input file

        Raises
        -------
        ValueError: if the input file has non-allowed extension
        """
        if self.check_extension():
            parsed_file = md.parse(self.path)
            logger.info(f'Successfully parsed file: {self.path}')
            return parsed_file
        else:
            logger.warning('The parsed file format is not valid, mzML or mzXML file is required.')
            raise argparse.ArgumentTypeError('Please parse an input file of either .mzML or .mzXML format.')

    def get_spectrum_list(self, parsed_file: xml.dom.minidom.Document) -> List[xml.dom.minidom.Element]:
        """Creates and returns list with spectra from the xml.dom.minidom.Document of the input .mzML file.

        Parameters
        ----------
        parsed_file: xml.dom.minidom.Document
            xml.dom.minidom.Document of the parsed input file

        Returns
        -------
        spectrum_list: List[xml.dom.minidom.Element]
            list containing the spectra as xml.dom.minidom.Element
        """
        if self.format == 'mzml':
            spectrum_list = parsed_file.getElementsByTagName('spectrum')
        elif self.format == 'mzxml':
            spectrum_list = parsed_file.getElementsByTagName('scan')
        if spectrum_list == list():
            logger.warning('The parsed file does not contain any spectrum data')
        return spectrum_list

    def get_spectrum_dict(self, spectrum_list: List[xml.dom.minidom.Element]) -> Dict[int, xml.dom.minidom.Element]:
        """Creates and returns dictionary with spectrum ids as key and xml.dom.minidom.Element as values from
        spectrum list.

        Parameters
        ----------
        spectrum_list: List[xml.dom.minidom.Element]
            list containing the spectra as xml.dom.minidom.Element

        Returns
        -------
        spectrum_dict: Dict[int, xml.dom.minidom.Element]
            dictionary with spectrum ids as key and xml.dom.minidom.Element as values

        """
        spectrum_dict = dict()
        if self.format == 'mzml':
            for spectrum in spectrum_list:
                ids = int(spectrum.getAttribute("index"))
                spectrum_dict[ids] = spectrum
        elif self.format == 'mzxml':
            for spectrum in spectrum_list:
                ids = int(spectrum.getAttribute("num"))
                spectrum_dict[ids] = spectrum
        return spectrum_dict

    def get_compression(self, spectrum_dict: Dict[int, xml.dom.minidom.Element]):
        """Gathers information about the encoding of the binary arrays in the parsed input file.

        Parameters
        ----------
        spectrum_dict: Dict[int, xml.dom.minidom.Element]
            dictionary with spectrum ids as key and xml.dom.minidom.Element as values

        Returns
        -------

        """
        if self.format == 'mzxml':
            logger.warning('Extracting compression types is only available for .mzML files.')
            raise argparse.ArgumentTypeError('Extracting compression types is only available for .mzML files.')
        compression_dict = dict()
        for key in spectrum_dict:
            params = spectrum_dict[key].getElementsByTagName('binaryDataArray')
            data = list()
            for array in params:
                p = array.getElementsByTagName('cvParam')
                for e in p:
                    data.append(e.getAttribute('name'))
            compression_dict[key] = data
        compression_list = list()
        for key in compression_dict:
            mz = compression_dict[key].index('m/z array')
            intensity = compression_dict[key].index('intensity array')
            length = len(compression_dict[key])
            if mz < intensity:
                mz_comp = compression_dict[key][:int(length/2)]
                int_comp = compression_dict[key][int(length/2):]
            else:
                int_comp = compression_dict[key][:int(length/2)]
                mz_comp = compression_dict[key][int(length/2):]
            compression_list.append([mz_comp, int_comp])
        dict_final = dict()
        for key in compression_dict:
            dict_final[key] = {'mz': dict(), 'intensity': dict()}
            value_mz = compression_list[key][0]
            value_int = compression_list[key][1]
            for val in value_mz:
                if 'float' in val or 'bit' in val:
                    dict_final[key]['mz']['data_type'] = val
                elif 'compression' in val:
                    dict_final[key]['mz']['compression'] = val
            for val in value_int:
                if 'float' in val or 'bit' in val:
                    dict_final[key]['intensity']['data_type'] = val
                elif 'compression' in val:
                    dict_final[key]['intensity']['compression'] = val
        self.compression = dict_final
        return

    def get_binary_spectrum_values(self, spectrum_dict: Dict[int, xml.dom.minidom.Element]):
        """Extracts binary data arrays for each spectrum.

        Parameters
        ----------
        spectrum_dict: Dict[int, xml.dom.minidom.Element]
            dictionary with spectrum ids as key and xml.dom.minidom.Element as values

        Returns
        -------

        """
        if self.format == 'mzxml':
            logger.warning('Binary data extraction is only available for .mzML files.')
            raise argparse.ArgumentTypeError('Binary data extraction is only available for .mzML files.')
        vals = dict()
        for key in spectrum_dict:
            if spectrum_dict[key].getAttribute('defaultArrayLength') != '0':
                mz_array, intensity_array = spectrum_dict[key].getElementsByTagName('binary')
                vals[key] = {'mz': mz_array.firstChild.nodeValue,
                             'intensity': intensity_array.firstChild.nodeValue}
            else:
                vals[key] = {'mz': None, 'intensity': None}
        self.binary_values = vals
        return

    def decode_decompress(self):
        """Takes the raw spectrum values and creates a dictionary of decoded and uncompressed m/z and intensity values."""
        if self.format == 'mzxml':
            logger.warning('Decoding binary arrays is only available for .mzML files.')
            raise argparse.ArgumentTypeError('Decoding binary arrays is only available for .mzML files.')
        spectrum_data = dict()
        for key in self.binary_values:
            encoded_mz_data, encoded_int_data = self.binary_values[key]['mz'], self.binary_values[key]['intensity']
            if encoded_mz_data is not None or encoded_int_data is not None:
                decoded_mz_data = base64.standard_b64decode(encoded_mz_data)
                decoded_int_data = base64.standard_b64decode(encoded_int_data)  # # decodes the string
                if self.compression[key]['mz']['compression'] == 'zlib compression':
                    decompressed_mz_data = zlib.decompress(decoded_mz_data)
                else:
                    decompressed_mz_data = decoded_mz_data
                if self.compression[key]['intensity']['compression'] == 'zlib compression':
                    decompressed_int_data = zlib.decompress(decoded_int_data)  # decompresses the data
                else:
                    decompressed_int_data = decoded_int_data
                if self.compression[key]['mz']['data_type'] == '32-bit float':
                    mz_data = struct.unpack('<%sf' % (len(decompressed_mz_data) // 4),
                                            decompressed_mz_data)  # unpacks 32-bit m/z values as floats
                elif self.compression[key]['mz']['data_type'] == '64-bit float':
                    mz_data = struct.unpack('<%sd' % (len(decompressed_mz_data) // 8),
                                            decompressed_mz_data)
                if self.compression[key]['intensity']['data_type'] == '32-bit float':
                    int_data = struct.unpack('<%sf' % (len(decompressed_int_data) // 4),
                                             decompressed_int_data)  # unpacks 32-bit intensity values as floats
                elif self.compression[key]['intensity']['data_type'] == '64-bit float':
                    int_data = struct.unpack('<%sd' % (len(decompressed_int_data) // 8),
                                             decompressed_int_data)
                spectrum_data[key] = {'mz': mz_data, 'intensity': int_data}
            else:
                spectrum_data[key] = {'mz': None, 'intensity': None}
        self.spectrum_data = spectrum_data
        return

    def get_values(self, spectrum_dict: Dict[int, xml.dom.minidom.Element]) -> Dict[int, Dict]:
        """Creates dictionary with spectrum ids and base peak m/z, base peak intensity, total ion current,
        lowest and highest observed m/z.

        Parameters
        ----------
        spectrum_dict: Dict[int, xml.dom.minidom.Element]
            dictionary with spectrum ids as key and xml.dom.minidom.Element as values

        Returns
        -------
        value_dict: Dict[int, Dict]
            dictionary containing spectrum ids and the spectrum values
        """
        value_dict = dict()
        if self.format == 'mzml':
            for index, key in enumerate(spectrum_dict):
                if spectrum_dict[key].getAttribute('defaultArrayLength') != '0':
                    vals = spectrum_dict[key].getElementsByTagName('userParam')
                    bpmz = vals[0].getAttribute('value')
                    bpi = vals[1].getAttribute('value')
                    tic = vals[2].getAttribute('value')
                    lomz = vals[3].getAttribute('value')
                    homz = vals[4].getAttribute('value')
                    value_dict[key] = {'spectra_id': int(index), 'base_peak_m/z': round(float(bpmz), 2),
                                       'base_peak_intensity': round(float(bpi), 2),
                                       'total_ion_current': round(float(tic), 2),
                                       'lowest_observed_m/z': round(float(lomz), 2),
                                       'highest_observed_m/z': round(float(homz), 2)}
                else:
                    value_dict[key] = {'base_peak_m/z': None, 'base_peak_intensity': None, 'total_ion_current': None,
                                       'lowest_observed_m/z': None, 'highest_observed_m/z': None}
        elif self.format == 'mzxml':
            for index, key in enumerate(spectrum_dict):
                if spectrum_dict[key].getAttribute('peaksCount') >= '0':
                    bpmz = spectrum_dict[key].getAttribute('basePeakMz')
                    bpi = spectrum_dict[key].getAttribute('basePeakIntensity')
                    tic = spectrum_dict[key].getAttribute('totIonCurrent')
                    lomz = spectrum_dict[key].getAttribute('lowMz')
                    homz = spectrum_dict[key].getAttribute('highMz')
                    value_dict[key] = {'spectra_id': int(index), 'base_peak_m/z': round(float(bpmz), 2),
                                       'base_peak_intensity': round(float(bpi), 2),
                                       'total_ion_current': round(float(tic), 2),
                                       'lowest_observed_m/z': round(float(lomz), 2),
                                       'highest_observed_m/z': round(float(homz), 2)}
        self.values = value_dict
        logger.info('Successfully gathered spectrum values.')
        return value_dict

    def analyse_spectrum(self) -> pd.DataFrame:
        """Wrapper function for the parsing of an mzML file and the extraction of the m/z and intensity values.
        Returns
        -------
        df_values: pd.DataFrame
            Dataframe containing spectrum ids and base peak m/z, base peak intensity, total ion current,
            lowest and highest observed m/z.
        """
        parsed_file = self.parse_file()
        s_list = self.get_spectrum_list(parsed_file)
        spectrum_dictionary = self.get_spectrum_dict(s_list)
        if self.format == 'mzml':
            self.get_compression(spectrum_dictionary)
            self.get_binary_spectrum_values(spectrum_dictionary)
            self.decode_decompress()
        values_spectrum = self.get_values(spectrum_dictionary)
        df_values = pd.DataFrame.from_dict(values_spectrum, orient='index', columns=['spectra_id',
                                                                                     'base_peak_m/z',
                                                                                     'base_peak_intensity',
                                                                                     'total_ion_current',
                                                                                     'lowest_observed_m/z',
                                                                                     'highest_observed_m/z'])
        return df_values
