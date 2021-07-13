import requests
import shutil
import pandas as pd

from typing import Dict, List


def fetch_descartes_human_tissue(out_file: str, verbose: bool = True) -> None:
    """Function to fetch Loom Single-Cell tissue data from
    Descartes human database.

    Args:
        out_file: Output file that is going to store .loom data
        verbose: If True (default), print statements about download

    Examples:
        >>> fetch_descartes_human_tissue("Human_Tissue.loom")

    """
    url = (
        "https://shendure-web.gs.washington.edu/content/members/cao1025/"
        "public/FCA_RNA_supp_files/scanpy_cells_all/"
        "Human_RNA_processed.loom"
    )
    if verbose:
        print("Downloading Human Single-Cell data from Descartes database")
        print(f"data url: {url}")

    with requests.get(url, stream=True, timeout=60) as data:
        with open(out_file, 'wb') as out:
            shutil.copyfileobj(data.raw, out)
    if verbose:
        print(f"Downloaded data to {out_file}")


def fetch_descartes_by_tissue(
    list_tissues: List[str],
    out_dir: str,
    verbose: bool = True
) -> None:
    """Function to fetch Loom Single-Cell tissue data from
    Descartes human database by choosing which tissues will be donwloaded.

    Args:
        list_tissues: List of tissues names to be downloaded.
        out_dir: Output directory that is going to store .loom data.
        verbose: If True (default), print statements about download.

    Examples:
        >>> fetch_descartes_by_tissue(
                list_tissues=["Thymus", "Hearth"]
                out_dir="data"
            )

    """
    base_url = (
        "https://shendure-web.gs.washington.edu/content/members/cao1025/"
        "public/FCA_RNA_supp_files/scanpy_cells_by_tissue"
    )
    for tissue in list_tissues:
        url = f"{base_url}/{tissue}_processed.loom"
        if verbose:
            print((
                f"Downloading {tissue} tissue Human Single-Cell data "
                "from Descartes database"
            ))
            print(f"data url: {url}")

        file_name = f"{out_dir}/{tissue}_data.loom"
        with requests.get(url, stream=True, timeout=60) as data:
            with open(file_name, 'wb') as out:
                shutil.copyfileobj(data.raw, out)
        if verbose:
            print(f"Downloaded {file_name} to {out_dir}")


def fetch_de_genes_for_cell_type(
    verbose: bool = False
) -> Dict[str, List[str]]:
    """Function to fetch Differentially Expressed (DE) genes from Descartes
    Human Atlas from 77 Main Cell types found in 15 Organs.

    Args:
        verbose: If True (default), print statements about download

    Returns:
        Dictionary mapping each main cell type to its differentially
        expressed genes. Example: {
            "Acinar cells": ["MIR1302-11", "FAM138A", ...],
            "Myeloid cells": ["CU459201.1", "OR4G4P", ...] ...
        }

    """
    url = (
        "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/"
        "downloads/data_summarize_fetus_data/DE_gene_77_main_cell_type.csv"
    )
    if verbose:
        print((
            "Downloading Human Single-Cell Differentially Expressed"
            "genes for 77 Main Cell types found in 15 Organs."
        ))
        print(f"data url: {url}")

    de_df = pd.read_csv(url)

    cell_types = de_df["max.cluster"].unique()
    de_mapping = {}
    for type in cell_types:
        list_genes = de_df[
            de_df["max.cluster"] == type
        ]["gene_id"].tolist()
        list_genes = [gene.replace("'", "") for gene in list_genes]

        de_mapping[type] = list_genes

    return de_mapping
