import scanpy as sc
import pandas as pd
import upsetplot as upset

from anndata import AnnData
from typing import List
from matplotlib import pyplot as plt


def marker_genes(
    adata: AnnData,
    name: str = "marker_genes.pdf",
    out_dir: str = ".",
    plot_format: str = "dotplot",
    n_genes: int = 5
) -> None:
    """Plots marker genes found in scanpy rank_genes_groups function.

    Args:
        adata: AnnData structure with ranked genes for all the groups analyzed.
        name: Name of the plot file output.
        out_dir: Output directory to store plots.
        type: Type of plot, being possible all plots in
            scanpy.pl.rank_genes group.

    """
    plot_type = {
        "": sc.pl.rank_genes_groups,
        "violin": sc.pl.rank_genes_groups_violin,
        "stacked_volion": sc.pl.rank_genes_groups_stacked_violin,
        "heatmap": sc.pl.rank_genes_groups_heatmap,
        "dotplot": sc.pl.rank_genes_groups_dotplot,
        "matrixplot": sc.pl.rank_genes_groups_matrixplot,
        "tracksplot": sc.pl.rank_genes_groups_tracksplot
    }
    sc.settings.figdir = out_dir
    plot_type[plot_format](adata, n_genes=n_genes, save=name)


def pathways(adata: AnnData, cluster_name: str) -> pd.DataFrame:
    """Returns the pathways DataFrame from a cluster, creating a nice
    visualization tool of the pathways annotated in that cluster.

    Args:
        adata: AnnData structure with ranked genes for all the groups analyzed.
        cluster_name: Name of the cluster to be visualized.

    Returns:
        DataFrame, creating a nice visualization from it in Jupyter
        Notebook.

    """
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    return adata.uns["pathways"][cluster_name]


def shared_pathways(
    adata: AnnData,
    clusters: List[str] = [],
    file_name: str = "shared_pathways.png",
    dpi: int = 300,
    color: str = "cornflowerblue"
) -> pd.DataFrame:
    """Plot pathways shared between input clusters using UpSet.

    Args:
        adata: AnnData structure with ranked genes for all the groups analyzed.
        clusters: List of clusters names. Default: all clusters.
        file_name: Name of the plot file output. Default: shared_pathways
        dpi: Image DPI. Default: 300
        color: Matplotlib color of UpSet plot. Default: cornflowerblue

    """
    if not clusters:
        clusters = adata.uns["pathways"].keys()

    pathways = list(set(sum([
        list(adata.uns["pathways"][cluster_name].name)
        for cluster_name in clusters
    ], [])))

    presence_dict = {
        cluster_name: [
            1 if path_name in
            list(adata.uns["pathways"][cluster_name].name)
            else 0 for path_name in pathways
        ] for cluster_name in clusters
    }

    presence_dict["pathways"] = pathways

    presence_df = pd.DataFrame(presence_dict)

    names = list(presence_df.columns[:-1])
    presence = presence_df[names].astype(bool)
    presence = pd.concat(
        [
            presence,
            presence_df[
                [col for col in presence_df.columns if col not in presence]
            ]
        ],
        axis=1
    ).set_index(names)

    fig = plt.figure(dpi=dpi, figsize=(24, 16))
    upset.UpSet(
        presence,
        show_counts='%d',
        facecolor=color,
        sort_by='cardinality'
    ).plot(fig=fig)
    fig.savefig(file_name)
