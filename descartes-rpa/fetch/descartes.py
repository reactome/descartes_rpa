import requests


def fetch_descartes_tissue(tissue: str, out_file: str) -> None:
    """Function to fetch Raw Gene Count Sparse Matrices from Single-Cell
    Descartes human database.

    Args:
        tissue: Name of the tissue that is going to be downloaded
        out_file: Output file that is going to store .RDS data

    Examples:
        >>> fetch_descartes_tissue("Thymus", "Thymus.RDS")

    """
    url = (
        "https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/"
        f"data_summarize_fetus_data/{tissue}_gene_count.RDS"
    )
    data = requests.get(url, allow_redirects=True)
    with open(out_file, 'wb') as out:
        out.write(data.content)
