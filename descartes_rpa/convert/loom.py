import scanpy as sc

from anndata import AnnData


def loom_to_anndata(
    file: str,
) -> AnnData:
    """Function to transform loom file into AnnData Struct

    Args:
        file: .loom file name with tissue Single-Cell data from descartes

    Returns:
        Anndata object with Single-Cell data

    Examples:
        >>> loom_to_anndata(file="Thymus_processed.loom")

    """
    return sc.read_loom(file)
