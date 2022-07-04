# Group 3 (WS-2021)

## Programming Project 4 - Compute the Number of Peptides of a Given Total Mass

------------------

### Contributors

--------------

- Priya Kempanna
- Lara Schmalenstroer
- Rohitha Ravinder 
- Shubhi Ambast

## Introduction 

--------------

Mass spectrometry is used extensively in the biochemical field for calculating the molecular weight of isolated
proteins as well as identifying possible modifications made to the protein itself (post-translational modifications). MS works by
generating a spectrum (known as a mass spectrum) which plots the mass-to-charge ratio (m/z) of the ions in the sample. By
using the values derived from such a spectrum, one can determine with a high degree of accuracy what the compound (or
peptide) is.

![image](https://upload.wikimedia.org/wikipedia/commons/1/1f/Mass_spectrometry_protocol.png)
*https://upload.wikimedia.org/wikipedia/commons/1/1f/Mass_spectrometry_protocol.png*
## Repository skeleton

--------------
```
├── frontend
│ └── templates
    └── layout.html
    └── macros.html
    └── template.html
    └── upload.html
  └── run.py
├── mass_spectrum
│ └── ms_package
    └── __init__.py
    └── cli.py
    └── peptide_prediction.py
    └── protein_prediction.py
    └── reader.py
    └── startup.py
  └── tests
    └── data
    └── __init__.py
    └── constants.py
    └── test_peptide_prediction.py
    └── test_protein_prediction.py
    └── test_reader.py
  └── setup.py
├── Dockerfile
├── README.md
├── docker-compose.yml
├── Group03_contributions.pdf
├── project4_spectra.pdf
 
```

## Description about the package

---

ms_package includes python scripts and tests to perform following tasks:
1. **Parse raw MS files** for their m/z values using */reader.py*
2. **Predict which peptides** using */peptide_prediction.py*
3. **Determine the proteins** that the predicted peptides could be derived from using */protein_prediction.py*
4. **Create a frontend** to upload a raw MS output file and obtain a list of possible peptides
5. **Containerize the Application** to bundle backend and frontend together using Docker

- [pyOpenMS](https://pyopenms.readthedocs.io/en/latest/) is an open source Python library used in this project for analysis of mass spectrometry raw data(mzXML, mzML, TraML, fasta, pepxml). The package helps in identifying of peptide fragments, isotopic abundances and peptide search.
  - pyOpenMS interacts with other search engines such as Mascot, MSFragger, OMSSA, Sequest, SpectraST, XTandem to identify proteins from peptide sequence databases.

- [EMBL-EBI Proteins API](https://www.ebi.ac.uk/proteins/api/doc/#/proteomics) is used to make API query which searches in the databases like MaxQB, PeptideAtlas, EPD and  ProteomicsDB and saves the output in json format.
  - Base url API for the task used here: "https://www.ebi.ac.uk/proteins/api/proteomics?offset=0&size=100&peptide={api_query}" ; api_query = peptide(str)
  

## Dependencies

--------------

- [Numpy] (https://numpy.org/)
- [pyOpenMS] (https://pyopenms.readthedocs.io/en/latest/#)
- [Pandas] (https://pandas.pydata.org/)
- [Flask] (https://flask.palletsprojects.com/en/2.0.x/)
- [Click] (https://click.palletsprojects.com/en/8.0.x/)
- [Requests] (https://docs.python-requests.org/en/latest/)


## Installation

--------------

- To install pyOpenMS from command line:
```bash
pip install numpy
pip install pyopenms
```
- To install the package from command line, navigate to the folder package and execute:
```bash
pip install mass_spectrum
```
- Docker
```bash
- To install Docker on your machine - (https://docs.docker.com/get-docker/)
- Install the "Docker" plugin in PyCharm or whatever IDE you are using
```

## Usage

-------
```python
- import ms_package

    - ms_package get-spectrum-values /tests/data/BSA1.mzML -v

    - ms_package peptide-info /tests/data/BSA.fasta /tests/data/BSA1.mzML -v

    - ms_package protein-info -f /tests/data/BSA.fasta -m /tests/data/BSA1.mzML -v

```

```python
- to test dockerfile, navigate to the folder with the Dockerfile and execute :
    
    - $ docker build . -t ms:latest 
    - $ docker run --name ms_test -p 5000:5000 -d ms:latest

- to have Docker build the container execute the following command line from the directory with the docker-compose.yml file:
    - $ docker-compose up -d
```
## License

-------

Code released under the [MIT license] (https://choosealicense.com/licenses/mit/).

## References

-------

- [Accurate peak list extraction from proteomic mass spectra for identification and profiling studies](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-518)
- [Mascot Server](https://www.matrixscience.com/search_form_select.html)
- [Detecting protein variants by mass spectrometry: a comprehensive study in cancer cell-lines](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0454-9)
- [Protein identification and expression analysis using mass spectrometry](https://www.cell.com/trends/microbiology/fulltext/S0966-842X(06)00076-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0966842X0600076X%3Fshowall%3Dtrue)
