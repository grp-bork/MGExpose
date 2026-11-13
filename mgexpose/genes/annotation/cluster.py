# pylint: disable=R0913,R0914,R0917
""" Module to process gene clustering output. """

import logging

from contextlib import nullcontext

from .clustering_parser import (
    parse_db_clusters,
    parse_full_seq_clusters,
    parse_y_clusters,
    evaluate_y_clusters
)
from ..gene import Gene


logger = logging.getLogger(__name__)


def add_clusters_trec_format(fn, genes, genome_id=None, core_threshold=0.95,):
    """
    trec format cluster:
    (genome_id)_(gene_id), cluster_id, m=member_genomes n=pangenome_size m/n
    SAMEA112578415.psa_megahit.semibin.00344_k141_5555417_6 ANI_trec_47531.000096e  2       3       0.66667
    ------------------genome---------------- -----gene-----
    """

    n_genes, n_core_genes = 0, 0
    for n_genes, (gene_id, cluster_id, m, n, occ) in enumerate(parse_db_clusters(fn), start=1):

        gene = genes.get(gene_id.replace(genome_id + "_", ""))
        
        logger.info(
            "Checking cluster %s gene %s... %s",
            str(cluster_id), gene_id, str(gene),
        )
        if gene is not None:
            gene.is_core = float(occ) > core_threshold
            n_core_genes += gene.is_core
            gene.cluster = cluster_id

    logger.info(
        "Parsed %s precomputed gene-cluster mappings with %s core genes (%s%%)",
        n_genes,
        n_core_genes,
        f"{n_core_genes / n_genes:.3f}",
    )


def add_clusters(
    fn,
    genes,
    cluster_type=None,
    # use_y_clusters=False,
    core_threshold=0.95,
    output_dir=None,
    genome_id=None,
):
    """ Add information from gene clustering to allow for core/accessory gene classification """

    # ("y_spire", "y_denovo", "trec", "pg3_denovo", "pg3_precomputed",), default="y_denovo")
    # if use_y_clusters:
    if cluster_type[0] == "y":
        # if core_threshold == -1:
        if cluster_type == "y_denovo":
            parse_y_clusters(fn, genes)
        elif cluster_type == "y_spire":
            evaluate_y_clusters(fn, genes, core_threshold=core_threshold,)
        else:
            raise ValueError(f"Unknown {cluster_type=}")
        return None

    write_data = False
    gene_clusters_out = nullcontext()
    n_genomes = 0
    cluster_genes = {}

    with gene_clusters_out:
        if fn is not None:

            # if core_threshold != -1:
            if cluster_type == "pg3_denovo":
                n_genomes, cluster_genes, gene_cluster_map, _ = parse_full_seq_clusters(
                    genes,
                    fn,
                    output_dir=output_dir,
                )

                logger.info(
                    "Parsed %s genomes with %s gene clusters.",
                    n_genomes,
                    len(cluster_genes),
                )
            else:
                # "pg3_precomputed"
                gene_cluster_map = list(parse_db_clusters(fn))

                n_genes = len(gene_cluster_map)
                n_core_genes = sum(1 for _, _, is_core in gene_cluster_map if is_core)
                logger.info(
                    "Parsed %s precomputed gene-cluster mappings with %s core genes (%s%%)",
                    n_genes,
                    n_core_genes,
                    round(n_core_genes / n_genes, 2),
                )

            for gene_id, *cluster in gene_cluster_map:
                cluster, *is_core = cluster
                is_core = is_core[0].lower() == "true" if is_core else None
                if genome_id is None or gene_id.startswith(genome_id):
                    # i think this section is still influenced by old data formats?
                    # that's why i don't use genes.get_gene here for now!
                    if genome_id is None:
                        gene = genes.get(gene_id)
                    else:
                        gene = genes.get(gene_id.replace(genome_id + ".", ""))
                    # gene = genes.get(
                    #     gene_id,
                    #     genes.get(
                    #         gene_id.replace(genome_id + ".", "")
                    #     )
                    # )
                    # gene = genes.get_gene(gene_id)
                    logger.info(
                        "Checking cluster %s gene %s... %s",
                        str(cluster),
                        gene_id,
                        str(gene),
                    )
                    if gene is not None and gene.species is not None:
                        gene.cluster = cluster

                        if cluster_genes:
                            occ = cluster_genes[cluster]
                            gene.is_core = Gene.is_core_gene(
                                occ, n_genomes, core_threshold=core_threshold,
                            )
                        elif core_threshold == -1:
                            gene.is_core = is_core

                        if write_data:
                            print(gene.id, gene.cluster, sep="\t", file=gene_clusters_out)

    return None
