import scanpy as sc
import pandas as pd

from anndata import AnnData
from collections import OrderedDict
from typing import List
from reactome2py import analysis
from functools import reduce

from descartes_rpa.fetch.descartes import fetch_de_genes_for_cell_type


def list_to_str(gene_list: List[str]) -> str:
    """Transform list of genes into a large string separated by comma to
    feed reactome2py identifiers function

    Args:
        gene_list: List of genes short-names

    Returns:
        Large string with all genes short-names separated by comma

    """
    genes_str = ""
    for gene in gene_list:
        genes_str = genes_str + "," + gene
    return genes_str[1:]


def enrich_list(genes_str: str, species: str = "Homo Sapiens") -> pd.DataFrame:
    """Uses reactome2py identifiers to find enriched
    pathways associated with input list of genes and return it to further
    analysis

    Args:
        genes_str: Large string with all genes short-names separated by comma
        species: Species from the Single-Cell libraries.

    Returns:
        DataFrame with each pathway enriched found from input genes
        in Reactome

    """
    result = analysis.identifiers(
        ids=genes_str,
        interactors=False,
        page_size="1",
        page="1",
        species=species,
        sort_by="ENTITIES_FDR",
        order="ASC",
        resource="TOTAL",
        p_value="1",
        include_disease=False,
        min_entities=None,
        max_entities=None,
        projection=True
    )
    token = result["summary"]["token"]
    token_result = analysis.token(
        token,
        species=species,
        page_size="-1",
        page="-1",
        sort_by="ENTITIES_FDR",
        order="ASC",
        resource="TOTAL",
        p_value="1",
        include_disease=False,
        min_entities=None,
        max_entities=None
    )
    return pd.DataFrame(token_result["pathways"])


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


def get_pathways_for_group(
    adata: AnnData,
    groupby: str = "Main_cluster_name",
    species: str = "Homo Sapiens",
    method: str = "wilcoxon",
    n_genes: int = 25
) -> None:
    """Get enriched pathways from Reactome using reatome2py

    Args:
        adata: AnnData object from descartes atlas
        groupyby: Anndata.obs key to be used for gene grouping, such as
        'Main_cluster_name' or 'leiden'
        species: Species from the Single-Cell libraries.
        method: Methods available in scanpy.tl.rank_genes_group, such as
        't-test' or 'wilcoxon'
        n_genes: Number of genes to be used as markers for group enrichment

    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=n_genes
    )

    unique_groups = adata.obs[groupby].unique().to_list()
    pathway_dict = OrderedDict()
    for group in unique_groups:
        genes_str = list_to_str(
            gene_list=adata.uns["rank_genes_groups"]["names"]
            [group].tolist()
        )
        pathway_df = enrich_list(genes_str=genes_str)
        pathway_dict[group] = pathway_df

    adata.uns["pathways"] = pathway_dict


def enrich_de_cell_types(adata: AnnData) -> None:
    """Get enriched pathways from Differentially Expressed genes found in
    Main Cell types annotated by Descartes Human Atlas.

    Args:
        adata: AnnData object from descartes atlas

    """
    de_mapping = fetch_de_genes_for_cell_type()

    pathway_dict = OrderedDict()
    for cell_type, genes in de_mapping.items():
        de_genes = list_to_str(gene_list=genes)
        pathway_df = enrich_list(genes_str=de_genes)
        pathway_dict[cell_type] = pathway_df

    adata.uns["pathways_de"] = pathway_dict


def get_shared(adata: AnnData, clusters: List[str]) -> pd.DataFrame:
    """Returns pathways IDs and names shared between input clusters.

    Args:
        adata: AnnData structure with ranked genes for all the groups analyzed.
        clusters: List of clusters names.

    Returns:
        DataFrame with pathways shared between input clusters.

    """
    clusters_df = [adata.uns["pathways"][name] for name in clusters]
    shared_df = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on=['stId'],
            how='inner',
            suffixes=("", "_b")
        ),
        clusters_df
    )

    return shared_df[["stId", "name"]]
