import click
import ms_package.startup
from ms_package.reader import Reader
from ms_package.peptide_prediction import PeptideSearch
from ms_package.protein_prediction import ProteinSearch
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@click.group()
def main():
    """Entry method"""
    pass


@main.command()
@click.argument('path')
@click.option('-v', '--verbose', default=False, is_flag=True, help="When used, will print the paths to STDOUT.")
def get_spectrum_values(path: str, verbose: bool = False):
    """Generates dataframe consisting of the spectrum values from the input mzml/mzxml file."""
    reader = Reader(path=path)
    data = reader.analyse_spectrum()
    if verbose:
        click.echo(data)


@main.command()
@click.argument('fasta_path')
@click.argument('mzml_path')
@click.option('-v', '--verbose', default=False, is_flag=True, help='When used, will print to STDOUT.')
def peptide_info(fasta_path: str, mzml_path: str, verbose: bool = False):
    """Generates dataframe consisting of peptide properties and list of peptide hit sequences"""
    search = PeptideSearch(fasta_path=fasta_path, mzml_path=mzml_path)
    info = search.peptide_wrapper()[0]
    if verbose:
        click.echo(info)


@main.command()
@click.option('-f', '--fasta', default=None, help='FASTA file of protein, submitted along MZML file')
@click.option('-m', '--mzml', default=None, help='MZML file containing spectrum information, submitted along FASTA file')
@click.option('-p', '--peptide', default=None, help='List of peptides to map to proteins')
@click.option('-v', '--verbose', default=False, is_flag=True, help='When used, prints table to STDOUT.')
@click.option('-s', '--sequence', default=False, is_flag=True, help='Option to print protein sequence.')
@click.option('-o', '--output', default=None, help='File path to save protein information')
def protein_info(fasta: str, mzml: str, peptide: list,  output: str, verbose: bool = False, sequence: bool = False):
    """Generates dataframe of peptide mapping to get proteins.
    """
    info = None
    ans_with_seq = None
    ans_without_seq = None

    if peptide and (fasta and mzml):
        raise ImportError("Please load either a Peptide list OR a FASTA and MZML file, not all 3!")

    if fasta and mzml:
        pep_search = PeptideSearch(fasta_path=fasta, mzml_path=mzml)
        info = pep_search.peptide_wrapper()[0]
        peptide_list = pep_search.peptide_wrapper()[1]
        pro_search = ProteinSearch(peptide_list)
        pro_search.get_proteins()
        ans_with_seq = pro_search.ans_df
        ans_without_seq = ans_with_seq.drop('Sequence', axis=1)

    if peptide:
        pro_search = ProteinSearch(peptide)
        pro_search.get_proteins()
        ans_with_seq = pro_search.ans_df
        ans_without_seq = pro_search.ans_df.drop('Sequence', axis=1)

    if verbose:
        if info is not None:
            click.echo(info)
        if sequence:
            click.echo(ans_with_seq)
        elif not sequence:
            click.echo(ans_without_seq)

    if ans_with_seq is not None and ans_without_seq is not None:
        if output:
            if sequence:
                ans_with_seq.to_csv(output, index=False)
            elif not sequence:
                ans_without_seq.to_csv(output, index=False)


if __name__ == '__main__':
    main()
