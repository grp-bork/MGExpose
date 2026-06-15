""" Module to annotate gene sets. """

import os

from functools import partial

from .cluster import add_clusters
from .eggnog import add_eggnog_annotation, parse_emapper
from ..geneset import GeneSet
from .phage import PhageDetection
from .recombinases import add_recombinases
from .conjugation import add_conjugation_systems

from ...utils.readers import (
    parse_macsyfinder_report,
    read_recombinase_hits,
)


def annotate_genes(args):

    genes = GeneSet.from_file(
        args.input_genes,
        genome_id=args.genome_id,
        species=args.species,
        gene_type=args.input_gene_type,
        composite_gene_ids=args.dbformat != "PG3",
    )

    annotations, has_clusters = compile_annotations(args, multi_run=False,)

    gene_info_out = open(
        os.path.join(
            args.output_dir,
            f"{args.genome_id}.gene_info.txt",
        ),
        "wt",
        encoding="UTF-8",
    )

    gene_info_gff = open(
        os.path.join(
            args.output_dir,
            f"{args.genome_id}.gene_info.gff3",
        ),
        "wt",
        encoding="UTF-8",
    )

    with gene_info_out, gene_info_gff:
        annotated_genes = genes.annotate(
            annotations, stream=gene_info_out, gffstream=gene_info_gff,
        )

    return annotated_genes, has_clusters


def compile_annotations(args, multi_run=False,):
    """ Compile annotation functions according to input parameters. """
    annotations = []
    has_clusters = False

    if hasattr(args, "recombinases") and args.recombinases:
        recombinases = read_recombinase_hits(args.recombinases,)
        if multi_run:
            recombinases = list(recombinases)

        annotations.append(
            partial(
                add_recombinases,
                recombinases=recombinases,
            )
        )

    if hasattr(args, "conjugation_data") and args.conjugation_data:
        conjugation_systems = parse_macsyfinder_report(
            args.conjugation_data, args.conjugation_rules,
        )
        if multi_run:
            conjugation_systems = list(conjugation_systems)

        annotations.append(
            partial(
                add_conjugation_systems,
                conjugation_systems=conjugation_systems,
            )
        )

    if hasattr(args, "phage_and_cargo_data") and args.phage_and_cargo_data:
        eggnog_annotations = parse_emapper(
            args.phage_and_cargo_data,
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

    if hasattr(args, "cluster_data") and args.cluster_data:
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
        has_clusters = True

    # try:
    #     if args.recombinase_hits:

    #         recombinases = read_recombinase_hits(args.recombinases,)
    #         if multi_run:
    #             recombinases = list(recombinases)

    #         annotations.append(
    #             partial(
    #                 add_recombinases,
    #                 recombinases=recombinases,
    #             )
    #         )
    # except AttributeError as err:
    #     print(f"ERR: {err}")

    # try:
    #     if args.conjugation_data:

    #         conjugation_systems = parse_macsyfinder_report(
    #             args.conjugation_data, args.conjugation_rules,
    #         )
    #         if multi_run:
    #             conjugation_systems = list(conjugation_systems)

    #         annotations.append(
    #             partial(
    #                 add_conjugation_systems,
    #                 conjugation_systems=conjugation_systems,
    #             )
    #         )
    # except AttributeError as err:
    #     print(f"ERR: {err}")

    # try:
    #     if args.phage_eggnog_data:

    #         eggnog_annotations = parse_emapper(
    #             args.phage_eggnog_data,
    #             phage_annotation=PhageDetection(args.phage_filter_terms),
    #         )
    #         if multi_run:
    #             eggnog_annotations = list(eggnog_annotations)

    #         annotations.append(
    #             partial(
    #                 add_eggnog_annotation,
    #                 eggnog_annotations,
    #             )
    #         )
    # except AttributeError as err:
    #     print(f"ERR: {err}")

    # try:
    #     if args.cluster_data:
    #         annotations.append(
    #             partial(
    #                 add_clusters,
    #                 args.cluster_data,
    #                 use_y_clusters=('use_y_clusters' in args and args.use_y_clusters),
    #                 core_threshold=('core_threshold' in args and args.core_threshold) or 0.95,
    #                 output_dir=args.output_dir,
    #                 genome_id=args.genome_id,
    #             )
    #         )
    #     has_clusters = True
    # except AttributeError as err:
    #     print(f"ERR: {err}")

    print(f"{annotations=} {has_clusters=}")

    return annotations, has_clusters
