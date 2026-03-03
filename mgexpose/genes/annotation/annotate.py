""" Module to annotate gene sets. """

import os

from functools import partial

from .cluster import add_clusters
from .eggnog import add_eggnog_annotation
from .recombinases import add_recombinases
from .secretion import add_secretion_systems


def compile_annotations(args):
    """ Compile annotation functions according to input parameters. """
    annotations = []
    try:
        if args.recombinase_hits:
            annotations.append(
                partial(
                    add_recombinases,
                    fn=args.recombinase_hits,
                    pyhmmer=not args.pyhmmer,
                )
            )
    except AttributeError:
        pass

    try:
        if args.txs_macsy_report and args.txs_macsy_rules:
            annotations.append(
                partial(
                    add_secretion_systems,
                    args.txs_macsy_report,
                    args.txs_macsy_rules,
                )
            )
    except AttributeError:
        pass

    try:
        if args.phage_eggnog_data and args.phage_filter_terms:
            annotations.append(
                partial(
                    add_eggnog_annotation,
                    args.phage_eggnog_data,
                    args.phage_filter_terms,
                )
            )
    except AttributeError:
        pass

    try:
        if args.cluster_data:
            annotations.append(
                partial(
                    add_clusters,
                    args.cluster_data,
                    use_y_clusters=('use_y_clusters' in args and args.use_y_clusters),
                    core_threshold=('core_threshold' in args and args.core_threshold) or 0.95,
                    output_dir=args.output_dir,
                    genome_id=args.genome_id,
                )
            )
    except AttributeError:
        pass

    return annotations


def annotate_genes(genes, args,):
    """ Annotate genes with MGE-relevant data. """

    print("ARGS", args)

    for f_ann in compile_annotations(args):
        f_ann(genes)

    with open(
        os.path.join(
            args.output_dir,
            f"{args.genome_id}.gene_info.txt",
        ),
        "wt",
        encoding="UTF-8",
    ) as gene_info_out:
        genes.dump(gene_info_out)

    return list(genes.values())
