import requests
import asyncio
from typing import List
from aiohttp import ClientSession


def ensembl_to_gene(ensembl_ids: List[str]) -> List[str]:
    """Wrapper to async transform list of Ensembl IDs
    into list of HGNC gene names.

    Args:
        ensembl_ids: List of Ensembl IDs

    Returns:
        List of gene names for each ID

    """
    loop = asyncio.get_event_loop()
    gene_name_list = loop.run_until_complete(
        retrieve_gene_name(ensembl_ids=ensembl_ids)
    )
    return gene_name_list


async def retrieve_gene_name(ensembl_ids: List[str]) -> List[str]:
    """Use Ensembl REST API to map Ensembl ID to HGNC gene names.

    Args:
        ensembl_ids: List of Ensembl IDs

    Returns:
        List of gene names for each ID

    """
    async with ClientSession() as _:
        gene_name_list = await asyncio.gather(
            *[retrieve_api(id) for id in ensembl_ids]
        )
    return gene_name_list


async def retrieve_api(ensembl_id: str) -> str:
    """Async wrapper to request Ensembl REST API, returning HGNC gene name.

    Args:
        ensembl_id: Ensembl ID

    Returns:
        Gene name of input Ensembl ID from HGNC database

    """
    url = (
        "https://rest.ensembl.org"
        f"/lookup/id/{ensembl_id}?"
        "external_db=HGNC;"
        "object_type=gene;"
        "species=homo_sapiens"
    )
    gene_request = requests.get(
        url,
        headers={"Content-Type": "application/json"}
    )

    if not gene_request.ok:
        gene_request.raise_for_status()
        return False

    return gene_request.json()["display_name"]
