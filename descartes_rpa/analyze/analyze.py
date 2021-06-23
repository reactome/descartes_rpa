import scanpy as sc

from anndata import AnnData
from collections import OrderedDict
from typing import List
from reactome2py import analysis


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


def enrich_list(genes_str: str) -> dict:
    """Uses reactome2py identifiers to find enriched
    pathways associated with input list of genes and return it to further
    analysis

    Args:
        genes_str: Large string with all genes short-names separated by comma

    Returns:
        Dictionary with each pathway enriched found from input genes
        in Reactome

    """
    result = analysis.identifiers(
        ids=genes_str,
        interactors=False,
        page_size="1",
        page="1",
        species="Homo Sapiens",
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
        species='Homo sapiens',
        page_size='-1',
        page='-1',
        sort_by='ENTITIES_FDR',
        order='ASC',
        resource='TOTAL',
        p_value='1',
        include_disease=False,
        min_entities=None,
        max_entities=None
    )
    return token_result


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
    method: str = "wilcoxon"
) -> None:
    """Get pathways from Reactome using reatome2py

    Args:
        adata: AnnData object from descartes atlas
        groupyby: Anndata.obs key to be used for gene grouping, such as
        'Main_cluster_name' or 'leiden'
        method: Methods available in scanpy.tl.rank_genes_group, such as
        't-test' or 'wilcoxon'

    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method
    )

    unique_groups = adata.obs[groupby].unique().to_list()
    pathway_dict = OrderedDict()
    for group in unique_groups:
        genes_str = list_to_str(
            gene_list=adata.uns["rank_genes_groups"]["names"]\
                [group].tolist()
        )
        pathway_data = enrich_list(genes_str=genes_str)
        pathway_dict[group] = {
            "pathway_data": pathway_data
        }
    adata.uns["pathways"] = pathway_dict
