import requests
import shutil


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
