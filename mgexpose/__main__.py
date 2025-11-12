#!/usr/bin/env python

# pylint: disable=R0912,R0914,R0915,R0913,R0917

""" Mobile genetic element annotation """

import contextlib
import gzip
import logging
import os
import pathlib

from .gene_annotator import GeneAnnotator
from .gene_calling import run_pyrodigal
from .handle_args import handle_args
from .island_processing import (
    generate_island_set,
    annotate_islands,
    evaluate_islands,
    prepare_precomputed_islands
)
from .islands import MgeGenomicIsland
from .readers.readers import read_fasta, read_prodigal_gff, read_mge_rules
from .gffio import read_genomic_islands_gff
from .recombinase_annotation import run_pyhmmer


MGE_TABLE_HEADERS = \
    ("is_tn",) + \
    MgeGenomicIsland.TABLE_HEADERS[1:6] + \
    MgeGenomicIsland.TABLE_HEADERS[8:14] + \
    ("mgeR", "name", "genes",)

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


def dump_islands(islands, out_prefix, db, write_genes=False, add_functional_annotation=False):
    """ dump genomic islands to intermediate gff """
    with open(
        f"{out_prefix}.unannotated_islands.gff3",
        "wt", encoding="UTF-8"
    ) as _out:
        print("##gff-version 3", file=_out)
        for island in sorted(islands, key=lambda isl: isl.contig):
            island.to_gff(
                _out, db, write_genes=write_genes,
                add_functional_annotation=add_functional_annotation,
                intermediate_dump=True,
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


def write_final_results(
    recombinase_islands,
    output_dir,
    genome_id,
    output_suffix,
    dbformat=None,
    write_tsv=True,
    write_gff=True,
    write_genes_to_gff=True,
    add_functional_annotation=False,
    genome_seqs=None,
):
    """ write final results """

    outstream = contextlib.nullcontext()
    gff_outstream = contextlib.nullcontext()

    out_prefix = os.path.join(
        output_dir,
        f"{genome_id}.{output_suffix}"
    )

    if write_tsv:
        outstream = open(
            f"{out_prefix}.txt",
            "wt",
            encoding="UTF-8",
        )
    if write_gff:
        gff_outstream = open(
            f"{out_prefix}.gff3",
            "wt",
            encoding="UTF-8",
        )

    # Sort the list of MGEGenomicIslands based on contig names
    sorted_islands = sorted(recombinase_islands, key=lambda isl: isl.contig)
    islands_by_contig = {}

    with outstream, gff_outstream:
        # TSV header
        if write_tsv:
            print(*MGE_TABLE_HEADERS, sep="\t", file=outstream)
        # GFF3 header
        if write_gff:
            print("##gff-version 3", file=gff_outstream)

        # Start recording the outputs
        for island in sorted_islands:
            islands_by_contig.setdefault(island.contig, []).append(island)
            # TSV: ignore gene-wise annotations; each line is recombinase island,
            # all gene IDs are stored in a gene_list column
            # assert genome_id == island.genome
            if write_tsv:
                island.to_tsv(outstream)
            # GFF3: add individual genes annotation;
            # parent lines are recombinase islands, children lines are genes
            # GFF3 parent term: recombinase island
            if write_gff:
                island.to_gff(
                    gff_outstream,
                    source_db=dbformat,
                    write_genes=write_genes_to_gff,
                    add_functional_annotation=add_functional_annotation,
                )

        if genome_seqs is not None:
            with gzip.open(
                f"{out_prefix}.ffn.gz", 
                "wt",
            ) as _out:
                for header, seq in read_fasta(genome_seqs):
                    seqid, *_ = header.split(" ")
                    for island in islands_by_contig.get(seqid, []):
                        attribs = island.get_attribs()
                        try:
                            del attribs["ID"]
                        except KeyError:
                            pass
                        try:
                            del attribs["name"]
                        except KeyError:
                            pass
                        attrib_str = ";".join(f"{item[0]}={item[1]}" for item in attribs.items() if item[1])
                        print(
                            f">{island.get_id()} {attrib_str}", seq[island.start - 1: island.end], sep="\n", file=_out
                        )

def annotate_genes(args, debug_dir=None,):
    annotator = GeneAnnotator(
        args.genome_id,
        args.speci,
        read_prodigal_gff(args.prodigal_gff),
        include_genome_id=args.include_genome_id,
        has_batch_data=args.allow_batch_data,
        dbformat=args.dbformat,
    )

    with open(
            os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt"),
            "wt",
            encoding="UTF-8",
    ) as gene_info_out:

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

        annotator.dump_genes(gene_info_out)

        return list(annotated_genes)



def denovo_annotation(args, debug_dir=None):
    """ denovo annotation """
    out_prefix = os.path.join(args.output_dir, args.genome_id)

    annotated_genes = annotate_genes(args, debug_dir=debug_dir,)
    # annotator = GeneAnnotator(
    #     args.genome_id,
    #     args.speci,
    #     read_prodigal_gff(args.prodigal_gff),
    #     include_genome_id=args.include_genome_id,
    #     has_batch_data=args.allow_batch_data,
    #     dbformat=args.dbformat,
    # )

    
    # with open(
    #         # os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt"),
    #         f"{out_prefix}.gene_info.txt",
    #         "wt",
    #         encoding="UTF-8",
    # ) as gene_info_out:        
    #     annotated_genes = annotator.annotate_genes(
    #         args.recombinase_hits,
    #         (
    #             args.phage_eggnog_data,
    #             args.phage_filter_terms,
    #         ),
    #         (
    #             args.txs_macsy_report,
    #             args.txs_macsy_rules,
    #             # args.macsy_version,
    #         ),
    #         clusters=args.cluster_data,
    #         use_y_clusters=args.use_y_clusters,
    #         core_threshold=(args.core_threshold, -1)[args.precomputed_core_genes],
    #         output_dir=args.output_dir,
    #         pyhmmer=args.pyhmmer_input,
    #         dump_gene_info=gene_info_out,
    #     )

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

        # 
        # test_islands = list(read_genomic_islands_gff(f"{out_prefix}.unannotated_islands.gff3"))
        # dump_islands(
        #     test_islands,
        #     out_prefix + ".test",
        #     args.dbformat,
        #     write_genes=True,
        #     add_functional_annotation=args.add_functional_annotation,
        # )

    # with open(
    #         os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt"),
    #         "wt",
    #         encoding="UTF-8",
    # ) as _out:
    #     annotator.dump_genes(_out)

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

    genomic_islands = None
    skip_island_identification = True

    if args.command == "denovo":
        skip_island_identification = args.skip_island_identification
        genomic_islands = denovo_annotation(args, debug_dir=debug_dir)

    elif args.command == "annotate_genes":
        annotate_genes(args, debug_dir=debug_dir,)

    elif args.command == "call_genes":
        run_pyrodigal(args)

    elif args.command == "annotate_recombinases":
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
