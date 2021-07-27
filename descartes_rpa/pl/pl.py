import scanpy as sc
import pandas as pd

from anndata import AnnData


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
