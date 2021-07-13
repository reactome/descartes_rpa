from anndata import AnnData


def save_data_with_pathways(
    adata: AnnData,
    directory: str,
    file: str,
    verbose: bool = True
) -> None:
    """Saves AnnData structure data and pathways annotated metadata.

    Args:
        adata: AnnData object from descartes atlas
        directory: Directory with pathways metadata annotated and
            AnnData struct
        file: Filename of h5ad AnnData file.
        verbose: If True (default), print statements about saving files

    Examples:
        >>> save_data_with_pathways(adata, directory="output", file="Thymus")

    """
    pathways_data = adata.uns["pathways"].copy()
    del adata.uns["pathways"]

    if verbose:
        print(f"Saving AnnData structure to {directory}/{file}.h5ad")

    adata.write(f"{directory}/{file}.h5ad")

    for pathway, df in pathways_data.items():
        pathway_name = pathway.replace(" ", "_")
        if verbose:
            print(
                f"Saving pathway data from {pathway} clusters to {directory}"
            )
        df.to_csv(f"{directory}/{pathway_name}_pathways.csv", index=False)

    adata.uns["pathways"] = pathways_data
