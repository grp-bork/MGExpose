""" Module for argument handling """

import argparse

from .readers.eggnog import EggnogReader


from . import __version__

def handle_args():
    """ Argument handling """
    ap = argparse.ArgumentParser(
        prog="mgexpose",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    ap.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    # ap.add_argument("--output_dir", "-o", type=str, default=".")
    # ap.add_argument("--dbformat", type=str, choices=("PG3", "SPIRE"))
    # ap.add_argument("--write_gff", action="store_true")
    # ap.add_argument("--write_genes_to_gff", action="store_true")
    # ap.add_argument("--dump_intermediate_steps", action="store_true")
    # ap.add_argument("--output_suffix", type=str, default="full_length_MGE_assignments")
    # ap.add_argument("--debug", action="store_true")

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
    denovo_ap.add_argument("prodigal_gff", type=str)
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

    return ap.parse_args()


def handle_args_old():
    """ Argument handling """
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_id", type=str)
    ap.add_argument("prodigal_gff", type=str)
    ap.add_argument("recombinase_hits", type=str)
    ap.add_argument("speci", type=str)

    ap.add_argument(
        "txs_macsy_rules",
        type=str,
        help=(
            "In macsyfinder v1, this is found in macsyfinder.summary(.txt)."
            " In v2+, this is provided with the pipeline."
        ),
    )
    ap.add_argument("txs_macsy_report", type=str)
    ap.add_argument("phage_eggnog_data", type=str)
    ap.add_argument("mge_rules", type=str)

    ap.add_argument("--cluster_data", type=str)
    ap.add_argument("--output_dir", "-o", type=str, default=".")
    ap.add_argument("--phage_filter_terms", type=str)
    ap.add_argument("--include_genome_id", action="store_true")
    ap.add_argument("--core_threshold", type=float, default=0.95)
    ap.add_argument("--macsy_version", type=int, choices=(1, 2), default=2)
    ap.add_argument(
        "--emapper_version",
        type=str,
        choices=EggnogReader.EMAPPER_FIELDS.keys(),
        default="v2.1.2",
    )
    ap.add_argument(
        "--allow_batch_data",
        action="store_true",
        help=(
            "SPIRE annotation may have data that does not relate to the current bin."
            " Ignore those data."
        ),
    )
    ap.add_argument(
        "--use_y_clusters",
        action="store_true",
        help=(
            "Gene clustering is performed against annotated"
            " and redundancy-reduced reference sets."
        ),
    )
    ap.add_argument(
        "--single_island",
        action="store_true",
        help="Input is genomic region, skips island computation."
    )
    ap.add_argument(
        "--precomputed_islands",
        type=str,
        help="Input is set of genomic regions, skips island computation."
    )
    ap.add_argument("--write_gff", action="store_true")
    ap.add_argument("--write_genes_to_gff", action="store_true")
    ap.add_argument("--add_functional_annotation",
                    action="store_true",
                    help="If specified, per gene emapper annotations are stored in the gff.")
    # ensure newest eggnog version
    ap.add_argument("--dump_intermediate_steps", action="store_true")
    ap.add_argument("--output_suffix", type=str, default="full_length_MGE_assignments")
    ap.add_argument("--dbformat", type=str, choices=("PG3", "SPIRE"))
    ap.add_argument(
        "--precomputed_core_genes",
        action="store_true",
        help="Core/accessory gene sets were precomputed."
    )
    ap.add_argument("--skip_island_identification", action="store_true")
    ap.add_argument("--extract_islands", type=str)

    return ap.parse_args()
