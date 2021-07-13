import pandas as pd
import scanpy as sc

from anndata import AnnData
from collections import OrderedDict
from glob import glob


def load_data_with_pathways(
    directory: str,
    verbose: bool = True
) -> AnnData:
    """Loads AnnData structure data being in expected directory format.

    Args:
        directory: Directory with pathways metadata annotated and
            AnnData struct
        verbose: If True (default), print statements about saving files

    Returns:
        AnnData with descartes data and pathway annotations from Reactome

    Examples:
        >>> adata = load_data_with_pathways(directory="output")

    """
    files = glob(f"{directory}/*")

    h5ad_file = [file for file in files if "h5ad" in file][0]

    pathways = OrderedDict()
    clusters = [file for file in files if ".csv" in file]
    for cluster in clusters:
        if verbose:
            print(f"Loading {cluster} pathway data.")
        cluster_name = cluster.split("_pathways")[0].split("/")[-1]
        pathways[cluster_name] = pd.read_csv(cluster)

    print(f"Loading {h5ad_file} AnnData file.")
    adata = sc.read(h5ad_file)
    adata.uns["pathways"] = pathways
    return adata
