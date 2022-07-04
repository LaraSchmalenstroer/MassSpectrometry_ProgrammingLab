import os
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Default Directory Paths
# home_dir = os.path.expanduser("~")  # Also acceptable
home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".plab2_group3_project")
LOG_DIR = PROJECT_DIR.joinpath("logs")
DATA_DIR = PROJECT_DIR.joinpath("data")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)


# Logging Configuration
LOG_FILE_PATH = os.path.join(LOG_DIR, "group3's_log.log")
logging.basicConfig(filename=LOG_FILE_PATH, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

