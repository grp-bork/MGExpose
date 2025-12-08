""" Module for argument handling """

import argparse
import logging

from .eggnog import EggnogReader

from . import __version__


def handle_args(args):
    """ Argument handling """

    log_ap = argparse.ArgumentParser(prog="mgexpose", add_help=False)
    log_ap.add_argument("-l", "--log_level", type=int, choices=range(6), default=logging.INFO)
    log_args, _ = log_ap.parse_known_args(args)

    try:
        logging.basicConfig(
            level=10 * log_args.log_level,
            format='[%(asctime)s] %(message)s'
        )
    except ValueError as invalid_loglevel_err:
        raise ValueError(f"Invalid log level: {log_args.log_level}") from invalid_loglevel_err

    ap = argparse.ArgumentParser(
        prog="mgexpose",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=(log_ap,),
    )

    ap.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    subparsers = ap.add_subparsers(dest="command", required=True)

    parent_subparser = argparse.ArgumentParser(add_help=False)
    parent_subparser.add_argument("--output_dir", "-o", type=str, default=".")
    parent_subparser.add_argument("--dbformat", type=str, choices=("PG3", "SPIRE"))
    parent_subparser.add_argument("--write_gff", action="store_true")
    parent_subparser.add_argument("--write_genes_to_gff", action="store_true")
    parent_subparser.add_argument("--dump_intermediate_steps", action="store_true")
    parent_subparser.add_argument(
        "--output_suffix", type=str, default="full_length_MGE_assignments",
    )
    parent_subparser.add_argument("--debug", action="store_true")

    denovo_ap = subparsers.add_parser(
        "denovo",
        help="Classify and annotate mobile genomic regions from annotated genes.",
        parents=(parent_subparser,),
    )
    denovo_ap.add_argument("genome_id", type=str)
    # denovo_ap.add_argument("prodigal_gff", type=str)
    denovo_ap.add_argument("input_genes", type=str)
    denovo_ap.add_argument("recombinase_hits", type=str)
    denovo_ap.add_argument("mge_rules", type=str)
    denovo_ap.add_argument("--speci", type=str, default="no_speci")
    denovo_ap.add_argument("--txs_macsy_rules", type=str)
    denovo_ap.add_argument("--txs_macsy_report", type=str)
    denovo_ap.add_argument("--phage_eggnog_data", type=str)
    denovo_ap.add_argument("--cluster_data", type=str)
    denovo_ap.add_argument("--skip_island_identification", action="store_true")
    denovo_ap.add_argument("--dump_genomic_islands", action="store_true")
    denovo_ap.add_argument("--phage_filter_terms", type=str)
    denovo_ap.add_argument("--input_gene_type", type=str, choices=("prodigal", "preannotated",), default="prodigal",)

    denovo_ap.add_argument("--include_genome_id", action="store_true")
    denovo_ap.add_argument("--core_threshold", type=float, default=0.95)
    denovo_ap.add_argument(
        "--allow_batch_data",
        action="store_true",
        help=(
            "SPIRE annotation may have data that does not relate to the current bin."
            " Ignore those data."
        ),
    )
    denovo_ap.add_argument(
        "--use_y_clusters",
        action="store_true",
        help=(
            "Gene clustering is performed against annotated"
            " and redundancy-reduced reference sets."
        ),
    )
    denovo_ap.add_argument(
        "--single_island",
        action="store_true",
        help="Input is genomic region, skips island computation."
    )
    denovo_ap.add_argument(
        "--precomputed_islands",
        type=str,
        help="Input is set of genomic regions, skips island computation."
    )
    denovo_ap.add_argument(
        "--precomputed_core_genes",
        action="store_true",
        help="Core/accessory gene sets were precomputed."
    )

    denovo_ap.add_argument(
        "--add_functional_annotation",
        action="store_true",
        help="If specified, per gene emapper annotations are stored in the gff."
    )
    # ensure newest eggnog version
    denovo_ap.add_argument("--extract_islands", type=str)

    denovo_ap.add_argument("--pyhmmer_input", action="store_true")

    denovo_ap.set_defaults(func=None)  # TODO

    identify_mobile_islands_ap = subparsers.add_parser(
        "identify_mobile_islands",
        help="Identify and classify genomic islands as mobile.",
        parents=(parent_subparser,),
    )

    identify_mobile_islands_ap.add_argument("island_gff", type=str)

    identify_mobile_islands_ap.set_defaults(func=None)  # TODO


    call_genes_ap = subparsers.add_parser(
        "call_genes",
        help="Call genes with Pyrodigal",
        parents=(parent_subparser,),
    )

    call_genes_ap.add_argument("genome_fasta", type=str)
    call_genes_ap.add_argument("genome_id", type=str)
    call_genes_ap.add_argument("--threads", "-t", type=int, default=1)
    call_genes_ap.set_defaults(func=None)  # TODO


    annotate_recombinases_ap = subparsers.add_parser(
        "recombinase_scan",
        help="Detect recombinases with PyHMMer",
        parents=(parent_subparser,),
    )

    annotate_recombinases_ap.add_argument("proteins_fasta", type=str)
    annotate_recombinases_ap.add_argument("gff", type=str)
    annotate_recombinases_ap.add_argument("recombinase_hmms", type=str)
    annotate_recombinases_ap.add_argument("mge_rules", type=str)
    annotate_recombinases_ap.add_argument("genome_id", type=str)
    annotate_recombinases_ap.add_argument("--threads", "-t", type=int, default=1)
    annotate_recombinases_ap.set_defaults(func=None)  # TODO

    return ap.parse_args()
