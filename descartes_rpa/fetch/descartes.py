import requests


def fetch_descartes_tissue(tissue: str, out_file: str) -> None:
    """Function to fetch Loom Single-Cell tissue data from
    Descartes human database.

    Args:
        tissue: Name of the tissue that is going to be downloaded
        out_file: Output file that is going to store .loom data

    Examples:
        >>> fetch_descartes_tissue("Thymus", "Thymus.loom")

    """
    url = (
        "https://shendure-web.gs.washington.edu/content/members/cao1025/"
        "public/FCA_RNA_supp_files/scanpy_cells_by_tissue/"
        f"{tissue}_processed.loom"
    )
    data = requests.get(url, allow_redirects=True)
    with open(out_file, 'wb') as out:
        out.write(data.content)
