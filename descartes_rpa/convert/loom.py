import scanpy as sc

from anndata import AnnData


def loom_to_anndata(
    file: str,
) -> AnnData:
    """Function to transform R .RDS file into a scipy sparse matrix

    Args:
        file: .loom file name with tissue Single-Cell data from descartes

    Examples:
        >>> loom_to_anndata(file="Thymus_processed.loom")

    """
    return sc.read_loom(file)