""" Module to annotate gene sets. """

from functools import partial

from .cluster import add_clusters
from .eggnog import add_eggnog_annotation, parse_emapper
from .phage import PhageDetection
from .recombinases import add_recombinases
from .secretion import add_secretion_systems

from ...utils.readers import (
    parse_macsyfinder_report,
    read_recombinase_hits,
)


def compile_annotations(args, multi_run=False,):
    """ Compile annotation functions according to input parameters. """
    annotations = []

    if getattr(args, "recombinase_hits"):

        recombinases = read_recombinase_hits(
            args.recombinase_hits, pyhmmer=args.pyhmmer_input,
        )
        if multi_run:
            recombinases = list(recombinases)

        annotations.append(
            partial(
                add_recombinases,
                recombinases=recombinases,
            )
        )
    
    if getattr(args, "txs_macsy_report") and getattr(args, "txs_macsy_rules"):

        secretion_systems = parse_macsyfinder_report(
            args.txs_macsy_report, args.txs_macsy_rules,
        )
        if multi_run:
            secretion_systems = list(secretion_systems)

        annotations.append(
            partial(
                add_secretion_systems,
                secretion_systems=secretion_systems,
            )
        )
    
    if getattr(args, "phage_eggnog_data") and getattr(args, "phage_filter_terms"):

        eggnog_annotations = parse_emapper(
            args.phage_eggnog_data,
            phage_annotation=PhageDetection(args.phage_filter_terms),
        )
        if multi_run:
            eggnog_annotations = list(eggnog_annotations)

        annotations.append(
            partial(
                add_eggnog_annotation,
                eggnog_annotations,
            )
        )

    if getattr(args, "cluster_data"):
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

    return annotations
