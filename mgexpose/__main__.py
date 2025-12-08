#!/usr/bin/env python

# pylint: disable=R0912,R0914,R0915,R0913,R0917

""" Mobile genetic element annotation """

import contextlib
import logging
import os
import pathlib
import sys

from .gene_annotator import GeneAnnotator
from .gene_calling import run_pyrodigal
from .gffio import read_prodigal_gff
from .handle_args import handle_args
from .island_processing import (
    generate_island_set,
    annotate_islands,
    evaluate_islands,
    prepare_precomputed_islands
)
from .readers import read_mge_rules, read_preannotated_genes
from .recombinase_scan import run_pyhmmer
from .writers import dump_islands, write_final_results


logger = logging.getLogger(__name__)


def process_islands(genes, genome_id, single_island=None, island_file=None, output_dir=None,):
    """ helper function to declutter main() """
    precomputed_islands = prepare_precomputed_islands(
        single_island=single_island,
        island_file=island_file,
        genome_id=genome_id,
    )

    if output_dir:
        pang_calls_out = open(
            os.path.join(
                output_dir,
                f"{genome_id}.pan_genome_calls.txt"),
            "wt",
            encoding="UTF-8",
        )

        islands_out = open(
            os.path.join(
                output_dir,
                f"{genome_id}.pan_genome_islands.txt",
            ),
            "wt",
            encoding="UTF-8",
        )

        raw_islands_out = open(
            os.path.join(
                output_dir,
                "..",  # temporary! this is only until i know if this is final output or not
                f"{genome_id}.pan_genome_islands_raw.txt",
            ),
            "wt",
            encoding="UTF-8",
        )
    else:
        pang_calls_out, islands_out, raw_islands_out = [contextlib.nullcontext() for _ in range(3)]

    with pang_calls_out, islands_out, raw_islands_out:
        yield from generate_island_set(
            genes,
            pang_calls_out=pang_calls_out,
            raw_islands_out=raw_islands_out,
            islands_out=islands_out,
            precomputed_islands=precomputed_islands,
        )


def identify_recombinase_islands(islands, genome_id, mge_rules, output_dir=None):
    """Identify MGE-islands according to a set of rules
     using various signals annotated in the corresponding gene set. """
    if output_dir:
        step1_out = open(
            os.path.join(
                output_dir,
                f"{genome_id}.assign_mge.step1.txt",
            ),
            "wt",
            encoding="UTF-8",
        )

        step2_out = open(
            os.path.join(
                output_dir,
                f"{genome_id}.assign_mge.step2.txt",
            ),
            "wt",
            encoding="UTF-8",
        )

        step3_out = open(
            os.path.join(
                output_dir,
                f"{genome_id}.assign_mge.step3.txt",
            ),
            "wt",
            encoding="UTF-8",
        )

    else:
        step1_out, step2_out, step3_out = [contextlib.nullcontext() for _ in range(3)]

    with step1_out:
        annotated_islands = list(annotate_islands(islands, outstream=step1_out))
    with step2_out, step3_out:
        return list(
            evaluate_islands(
                annotated_islands,
                read_mge_rules(mge_rules),
                outstream=step2_out,
                outstream2=step3_out
            )
        )


def annotate_genes(args, debug_dir=None,):

    if args.input_gene_type == "prodigal":
        genes = read_prodigal_gff(args.input_genes)
    else:
        genes = read_preannotated_genes(args.input_genes)

    annotator = GeneAnnotator(
        args.genome_id,
        args.speci,
        genes,
        include_genome_id=args.include_genome_id,
        has_batch_data=args.allow_batch_data,
        composite_gene_ids=args.dbformat != "PG3",
    )

    with open(
            os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt"),
            "wt",
            encoding="UTF-8",
    ) as gene_info_out:

        annotated_genes = annotator.annotate_genes(
            args.recombinase_hits if args.recombinase_hits else None,
            (
                args.phage_eggnog_data,
                args.phage_filter_terms,
            ) if args.phage_eggnog_data and args.phage_filter_terms else None,
            (
                args.txs_macsy_report,
                args.txs_macsy_rules,
            ) if args.txs_macsy_report else None,
            clusters=args.cluster_data if args.cluster_data else None,
            use_y_clusters=args.use_y_clusters,
            core_threshold=(args.core_threshold, -1)[args.precomputed_core_genes and not args.use_y_clusters],
            output_dir=args.output_dir,
            pyhmmer=args.pyhmmer_input,
        )

        annotator.dump_genes(gene_info_out)

        return list(annotated_genes)



def denovo_annotation(args, debug_dir=None):
    """ denovo annotation """
    out_prefix = os.path.join(args.output_dir, args.genome_id)

    annotated_genes = annotate_genes(args, debug_dir=debug_dir,)

    genomic_islands = list(
        process_islands(
            annotated_genes,
            args.genome_id,
            single_island=args.single_island,
            island_file=args.precomputed_islands,
            output_dir=debug_dir,
        )
    )

    if args.dump_genomic_islands or args.skip_island_identification:

        dump_islands(
            genomic_islands,
            out_prefix,
            args.dbformat,
            write_genes=True,
            add_functional_annotation=args.add_functional_annotation,
        )

    return genomic_islands


def main():
    """ main """

    args = handle_args(sys.argv[1:])
    logger.info("ARGS: %s", str(args))

    debug_dir = None
    cdir = args.output_dir
    if args.dump_intermediate_steps:
        cdir = debug_dir = os.path.join(args.output_dir, "debug")
    pathlib.Path(cdir).mkdir(exist_ok=True, parents=True)

    genomic_islands = None
    skip_island_identification = True

    if args.command == "denovo":
        skip_island_identification = args.skip_island_identification
        genomic_islands = denovo_annotation(args, debug_dir=debug_dir)

    elif args.command == "annotate_genes":
        annotate_genes(args, debug_dir=debug_dir,)

    elif args.command == "call_genes":
        run_pyrodigal(args)

    elif args.command == "recombinase_scan":
        run_pyhmmer(args)

    elif args.command == "annotate":
        raise NotImplementedError

    if not skip_island_identification:

        recombinase_islands = identify_recombinase_islands(
            genomic_islands,
            args.genome_id,
            args.mge_rules,
            output_dir=debug_dir,
        )

        if recombinase_islands:
            write_final_results(
                recombinase_islands,
                args.output_dir,
                args.genome_id,
                args.output_suffix,
                dbformat=args.dbformat,
                write_gff=args.write_gff,
                write_genes_to_gff=args.write_genes_to_gff,
                add_functional_annotation=args.add_functional_annotation,
                genome_seqs=args.extract_islands,
            )


if __name__ == "__main__":
    main()
