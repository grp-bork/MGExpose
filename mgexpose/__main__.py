#!/usr/bin/env python

# pylint: disable=R0912,R0914,R0915,R0913,R0917

""" Mobile genetic element annotation """

import contextlib
import logging
import os
import pathlib

from .gene_annotator import GeneAnnotator
from .gffio import write_island_gff
from .handle_args import handle_args
from .island_processing import (
    annotate_islands,
    evaluate_islands,
    generate_island_set,
    prepare_precomputed_islands
)
from .readers.readers import read_fasta, read_prodigal_gff, read_mge_rules
from .writers import write_final_results


logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s'
)

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

def denovo_annotation(args, debug_dir=None, contig_data=None,):
    """ denovo annotation """

    annotator = GeneAnnotator(
        args.genome_id,
        args.speci,
        read_prodigal_gff(args.prodigal_gff),
        include_genome_id=args.include_genome_id,
        has_batch_data=args.allow_batch_data,
        dbformat=args.dbformat,
    )

    annotated_genes = annotator.annotate_genes(
        args.recombinase_hits,
        (
            args.phage_eggnog_data,
            args.phage_filter_terms,
        ),
        (
            args.txs_macsy_report,
            args.txs_macsy_rules,
            # args.macsy_version,
        ),
        clusters=args.cluster_data,
        use_y_clusters=args.use_y_clusters,
        core_threshold=(args.core_threshold, -1)[args.precomputed_core_genes],
        output_dir=args.output_dir,
        pyhmmer=args.pyhmmer_input,
    )

    out_prefix = os.path.join(args.output_dir, args.genome_id)

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

        islands_by_contig = {}
        for island in sorted(genomic_islands, key=lambda isl: isl.contig):
            islands_by_contig.setdefault(island.contig, []).append(island)

        write_island_gff(
            f"{out_prefix}.unannotated_islands.gff3",
            islands_by_contig,
            contig_data,
            args.dbformat,
            write_genes_to_gff=True,
            add_functional_annotation=args.add_functional_annotation,
            intermediate_dump=True,
        )

        # dump_islands(
        #     genomic_islands,
        #     out_prefix,
        #     args.dbformat,
        #     write_genes=True,
        #     add_functional_annotation=args.add_functional_annotation,
        #     contig_data=contig_data,
        # )

        # 
        # test_islands = list(read_genomic_islands_gff(f"{out_prefix}.unannotated_islands.gff3"))
        # dump_islands(
        #     test_islands,
        #     out_prefix + ".test",
        #     args.dbformat,
        #     write_genes=True,
        #     add_functional_annotation=args.add_functional_annotation,
        # )

    with open(
            os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt"),
            "wt",
            encoding="UTF-8",
    ) as _out:
        annotator.dump_genes(_out)

    return genomic_islands


def main():
    """ main """

    args = handle_args()
    logger.info("ARGS: %s", str(args))

    debug_dir = None    
    cdir = args.output_dir
    if args.dump_intermediate_steps:
        cdir = debug_dir = os.path.join(args.output_dir, "debug")
    pathlib.Path(cdir).mkdir(exist_ok=True, parents=True)

    contig_data = {}
    if args.extract_islands:
        contig_data = {
            sid: len(seq)
            for sid, seq in read_fasta(args.extract_islands)
        }

    genomic_islands = None
    if args.command == "denovo":
        genomic_islands = denovo_annotation(args, debug_dir=debug_dir)

    elif args.command == "annotate":
        raise NotImplementedError

    if not args.skip_island_identification:

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
                contig_data=contig_data,
            )


if __name__ == "__main__":
    main()
