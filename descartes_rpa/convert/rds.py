import subprocess
import pandas as pd

from anndata import AnnData
from scipy.io import mmread
from scipy.sparse.coo import coo_matrix


def rds_to_anndata(
    file: str,
    r_script: str = "rds_to_mtx.R",
    matrix_file: str = "matrix.mtx",
    ensembl_ids: str = "ensembl_ids.csv",
    cell_ids: str = "cell_ids.csv"
) -> AnnData:
    """Function to transform R .RDS file into a scipy sparse matrix

    Args:
        file: .RDS file name with matrix data from descartes
        r_script: Rscript file that transforms .RDS to .MTX
        matrix_file: Output .MTX Sparse Matrix file
        ensembl_ids: .CSV file with ID for each gene in assay
        cell_ids: .CSV file with ID's for each cell in assay

    Examples:
        >>> rds_to_anndata(file="Thymus.RDS")

    """
    stdout, stderr = rds_to_mm(
        file=file,
        r_script=r_script,
        matrix_file=matrix_file,
        ensembl_ids=ensembl_ids,
        cell_ids=cell_ids
    )
    matrix = mmread(matrix_file)
    df = pd.DataFrame.sparse.from_spmatrix(matrix)
    genes = pd.read_csv(ensembl_ids)["x"]
    cells = pd.read_csv(cell_ids)["x"]
    df = df.set_index(genes)
    df.columns = cells
    df = df.T
    return AnnData(df)


def rds_to_mm(
    file: str,
    r_script: str,
    matrix_file: str,
    ensembl_ids: str,
    cell_ids: str
) -> coo_matrix:
    """Function to transform R .RDS file into a scipy sparse matrix

    Args:
        file: .RDS file with matrix data from descartes
        r_script: Rscript file that transforms .RDS to .MTX
        matrix_file: Output .MTX Sparse Matrix file
        ensembl_ids: .CSV file with ID for each gene in assay
        cell_ids: .CSV file with ID's for each cell in assay

    """
    r_cmd = (
        f"Rscript {r_script} "
        f"{file} "
        f"{matrix_file} "
        f"{ensembl_ids} "
        f"{cell_ids}"
    )
    process = subprocess.Popen(
        r_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True
    )
    stdout, stderr = process.communicate()
    stdout = stdout.decode("utf-8")
    stderr = stderr.decode("utf-8")
    return stdout, stderr
