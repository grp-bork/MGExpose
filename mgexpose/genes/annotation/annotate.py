""" Module to annotate gene sets. """

import os

from .cluster import add_clusters
from .eggnog import add_eggnog_annotation
from .recombinases import add_recombinases
from .secretion import add_secretion_systems

from ..geneset import GeneSet


def annotate_genes(args,):
    """ Annotate genes with MGE-relevant data. """

    genes = GeneSet.from_file(
        args.input_genes,
        genome_id=args.genome_id,
        speci=args.speci,
        gene_type=args.input_gene_type,
        composite_gene_ids=args.dbformat != "PG3",
    )

    if args.recombinase_hits is not None:
        add_recombinases(genes, args.recombinase_hits, pyhmmer=args.pyhmmer_input,)

    if args.txs_macsy_report is not None and args.txs_macy_rules is not None:
        add_secretion_systems(genes, args.txs_macsy_report, args.txs_macsy_rules,)

    if args.phage_eggnog_data is not None and args.phage_filter_terms is not None:
        add_eggnog_annotation(genes, args.phage_eggnog_data, args.phage_filter_terms,)

    if args.cluster_data is not None:
        add_clusters(
            genes,
            args.cluster_data,
            use_y_clusters=args.use_y_clusters,
            core_threshold=args.core_threshold,
            output_dir=args.output_dir,
        )

    with open(
        os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt"),
        "wt",
        encoding="UTF-8",
    ) as gene_info_out:
        genes.dump(gene_info_out)

    return list(genes.values())
