"""Test module constants."""
from pathlib import Path

TEST_FOLDER = Path(__file__).parent
TEST_DATA_DIR = TEST_FOLDER.joinpath("data")
TEST_MZML_FILE = TEST_DATA_DIR.joinpath("BSA1.mzML")
TEST_FASTA_FILE = TEST_DATA_DIR.joinpath("BSA.fasta")
TEST_MZXML_FILE = TEST_DATA_DIR.joinpath("7MIX_STD_110802_1.mzXML")

