import os
import sys

import logging
import requests

import json
import pandas as pd
from collections import defaultdict

from ms_package.startup import DATA_DIR

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ProteinSearch:
    """This class uses The Proteins API from EBI to map the given list of peptides."""

    def __init__(self, peptide_list: list):
        self.peptide_list = peptide_list

        self.sel_peptides = []
        self.protein_response = ""
        self.ans_df = None

    def filter_peptides(self) -> list:
        """Filters peptide sequence if it contains special characters.
        Returns
        ----------
        filtered_pep_list : list
            List of filtered peptide sequences.
        """
        filtered_pep_list = [pep for pep in self.peptide_list if "(" not in pep and ")" not in pep]
        return filtered_pep_list

    def divide_into_chunks(self, filtered_pep_list):
        """Divides the given list of peptide list into list of lists of 15 peptides
        to submit to API.
        Parameters
        ----------
        filtered_pep_list : list
            List of filtered peptide sequences.
        """

        logger.info(f"The peptide query consists of {len(filtered_pep_list)} peptides. "
                    f"If > 15, dividing it into chunks..")

        for i in range(0, len(filtered_pep_list), 15):
            chunks = filtered_pep_list[i:i + 15]
            self.sel_peptides.append(chunks)

    @staticmethod
    def get_api_query(list_peps: list) -> str:
        """Converts the list of peptides to The Proteins API query format.
        Parameters
        ----------
        list_peps: list
            List of peptides.

        Returns
        -------
        api_query : str
            Converted API format.
        """
        api_query = ""
        for pep in list_peps:
            if list_peps.index(pep) == 0:
                api_query = pep
            else:
                api_query += "%2C" + pep
        return api_query

    def proteins_api(self):
        """Make API query for 15 peptides at a time.
        The Proteins API is used, which searches in the databases: MaxQB, PeptideAtlas, EPD and  ProteomicsDB.
        The response in json format is saved.
        """

        logger.info("Calling The Proteins API...  "
                    "to match peptides in databases: MaxQB, PeptideAtlas, EPD and  ProteomicsDB")

        for pep_lil in self.sel_peptides:
            api_query = self.get_api_query(pep_lil)
            requestURL = f"https://www.ebi.ac.uk/proteins/api/proteomics?offset=0&size=100&peptide={api_query}"
            r = requests.get(requestURL, headers={"Accept": "application/json"})

            if not r.ok:
                r.raise_for_status()
                logger.error("API request not successful, Protein list not received")
                sys.exit()

            self.protein_response += r.text

        logger.info("Peptide mapping to Proteins completed successfully")

    def parse_content(self) -> pd.DataFrame:
        """Function to parse the API request into dataframe
        Returns:
        ----------
        dataframe :
            Dataframe of protein identification values.
        """

        ans_dict = defaultdict(lambda: defaultdict(str))
        for l in range(len(json.loads(self.protein_response))):
            ind = json.loads(self.protein_response)[l]

            ans_dict[l + 1]["Peptide"] = ind["features"][0]["peptide"]
            ans_dict[l + 1]["Protein_Accession"] = ind["accession"]
            ans_dict[l + 1]["Taxonomy_ID"] = ind["taxid"]

            source = []
            for data in ind['features'][0]["evidences"]:
                source.append(data["source"]["name"])

            ans_dict[l + 1]["Data Source"] = ",".join(source)
            ans_dict[l + 1]["Peptide Position"] = f"{ind['features'][0]['begin']} to {ind['features'][0]['end']}"

            ans_dict[l + 1]["Sequence"] = ind["sequence"]

        return pd.DataFrame.from_dict(ans_dict, orient="index")

    @staticmethod
    def file_handle():
        """Creates a directory and files to keep a track of peptide query submitted
        to the API.
        """
        pro_dir = os.path.join(DATA_DIR, "proteins")
        os.makedirs(pro_dir, exist_ok=True)

        pep_list_dir = os.path.join(pro_dir, "protein_data.csv")
        try:
            open(pep_list_dir, "x")
        except FileExistsError:
            pass

    def get_proteins(self):
        """ Wrapper function to call relevant methods.
        Only calls API if the data was not available from previous API queries.
        """
        self.file_handle()
        filtered = self.filter_peptides()

        pep_list_dir = os.path.join(DATA_DIR, "proteins", "protein_data.csv")

        try:
            df = pd.read_csv(pep_list_dir)
        except pd.errors.EmptyDataError:
            df = None

        if df is not None:
            exists_df = df[df['Peptide'].isin(filtered)]
            already_exists_peptide = list(exists_df["Peptide"])
            to_be_queried = [x for x in filtered if x not in already_exists_peptide]
        else:
            to_be_queried = filtered
            exists_df = None

        if to_be_queried:
            self.divide_into_chunks(to_be_queried)
            self.proteins_api()
            response_df = self.parse_content()
        else:
            response_df = None

        # half exists, half API
        if response_df is not None and exists_df is not None:
            self.ans_df = pd.concat([exists_df, response_df],  axis=0, ignore_index=True)
            response_df.to_csv(pep_list_dir, mode="a", header=False, index=False)

        # full API
        elif response_df is not None and exists_df is None:
            self.ans_df = response_df
            try:
                self.ans_df.to_csv(pep_list_dir, mode="a", index=False)
            except KeyError:
                self.ans_df.to_csv(pep_list_dir, mode="a", header=False, index=False)

        # full exists
        elif response_df is None and exists_df is not None:
            self.ans_df = exists_df
        else:
            logger.error("No response received")
