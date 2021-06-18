from anndata import AnnData


def scanpy_format(adata: AnnData) -> None:
    """From descartes loom data, format it to automatic plotting and API
    automatic analysis in scanpy

    Args:
        adata: AnnData object from descartes atlas

    """
    # Set gene short name as index
    adata.var.set_index("gene_short_name", inplace=True)
    adata.var_names_make_unique()

    # Change umap values to obsm, so scanpy plotting API can access
    adata.obsm["X_umap"] = adata.obs[
        ["Main_cluster_umap_1", "Main_cluster_umap_2"]
    ].to_numpy()
