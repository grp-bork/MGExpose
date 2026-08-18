#!/usr/bin/env python

# pylint: disable=R0912,R0914,R0915,R0913,R0917

""" Mobile genetic element annotation """

import logging
import os
import pathlib
import sys

from .genes.annotation.annotate import compile_annotations
from .genes.geneset import GeneSet
from .handle_args import handle_args
from .islands.annotated_genomic_island import AnnotatedGenomicIsland
from .islands.genomic_island import GenomicIsland
from .islands.island_processing import annotate_islands, evaluate_islands, prepare_islands
from .islands.mge_genomic_island import MgeGenomicIsland
from .utils.gffio import read_genomic_islands_gff, read_mge_genomic_islands_gff
from .utils.readers import read_precomputed_islands, read_mge_rules
from .utils.writers import extract_mge_seqs, write_final_results
from .utils.upstream.gene_calling import run_pyrodigal
from .utils.upstream.recombinase_scan import run_pyhmmer


logger = logging.getLogger(__name__)


def main():
    """ main """

    args = handle_args(sys.argv[1:])
    logger.info("ARGS: %s", str(args))

    cdir = args.output_dir
    if args.dump_intermediate_steps:
        cdir = os.path.join(args.output_dir, "debug")
    pathlib.Path(cdir).mkdir(exist_ok=True, parents=True)

    if args.command in ("call_genes", "recombinase_scan",):

        if args.command == "call_genes":
            run_pyrodigal(args)

        elif args.command == "recombinase_scan":
            run_pyhmmer(args)

    else:

        genomic_islands = None
        skip_island_identification = True

        if args.command == "denovo":
            skip_island_identification = args.skip_island_identification
            genes = GeneSet.from_file(
                args.input_genes,
                genome_id=args.genome_id,
                speci=args.speci,
                gene_type=args.input_gene_type,
                composite_gene_ids=args.dbformat != "PG3",
            )

            annotations = compile_annotations(args, multi_run=False,)

            with open(
                os.path.join(
                    args.output_dir,
                    f"{args.genome_id}.gene_info.txt",
                ),
                "wt",
                encoding="UTF-8",
            ) as gene_info_out:

                annotated_genes = genes.annotate(annotations, stream=gene_info_out,)

                precomputed_islands = None
                if args.single_island or args.precomputed_islands:
                    precomputed_islands = read_precomputed_islands(
                        genome_id=args.genome_id,
                        single_island=args.single_island,
                        island_file=args.precomputed_islands,
                    )

                genomic_islands = prepare_islands(
                    annotated_genes, precomputed_islands=precomputed_islands, contigs_are_islands=args.contigs_are_islands,
                )

                annotated_islands = annotate_islands(genomic_islands,)
                mge_islands = list(
                    evaluate_islands(annotated_islands, read_mge_rules(args.mge_rules),)
                )

                if mge_islands:
                    write_final_results(mge_islands, args)

        elif args.command == "reannotate":

            if args.annotation_mode == "mges":
                genomic_islands = read_mge_genomic_islands_gff(args.input_gff)
            elif args.annotation_mode == "raw_islands":
                genomic_islands = read_genomic_islands_gff(args.input_gff)
            # genomic_islands = read_island_gff(args.input_gff, GenomicIsland)
            mge_rules = read_mge_rules(args.mge_rules)
            out_prefix = os.path.join(
                args.output_dir,
                f"{args.genome_id}.mge_islands"
            )

            gene_info_out = open(
                os.path.join(args.output_dir, f"{args.genome_id}.gene_info.txt",),
                "wt", encoding="UTF-8",
            )

            gff_out = open(f"{out_prefix}.gff3", "wt", encoding="UTF-8",)

            with gene_info_out, gff_out:
                print("##gff-version 3", file=gff_out)
                mge_islands = {}

                annotations = compile_annotations(args, multi_run=True,)

                for ct, island in enumerate(genomic_islands):
                    # strip
                    island = GenomicIsland.from_island(island, genome_id=args.genome_id,)
                    genes = GeneSet.from_island(
                        island, composite_gene_ids=args.annotation_mode == "raw_islands",
                    )
                    # annotated_genes = annotate_genes(
                    # genes, annotations, gene_info_out, gene_info_header,)
                    annotated_genes = genes.annotate(
                        annotations, gene_info_out, with_header=ct == 0,
                    )
                    island.update_recombinases()
                    # annotated_island = AnnotatedGenomicIsland.from_genomic_island(island)
                    annotated_island = AnnotatedGenomicIsland.from_island(island)
                    # mge_island = MgeGenomicIsland.from_annotated_genomic_island(annotated_island)
                    mge_island = MgeGenomicIsland.from_island(annotated_island)
                    # mge_island = island
                    mge_island.evaluate_recombinases(mge_rules)
                    mge_island.to_gff(
                        gff_out,
                        source_db=None,
                    )
                    # mge_islands.append(mge_island)
                    mge_islands.setdefault(mge_island.contig, []).append(mge_island)
                if args.extract_islands:
                    extract_mge_seqs(args.extract_islands, mge_islands, out_prefix)

        elif args.command == "liftover":
            mge_islands = {}
            mge_rules = read_mge_rules(args.mge_rules)
            # with open(args.island_mapping, "rt") as _in:
            #     # island mapping format:
            #     # source -> dest, hence we need to reverse
            #     # as we're checking the dest islands only
            #     island_mapping = dict(
            #         line.strip().split()[::-1]
            #         for line in _in
            #     )
            island_mapping = dict([args.island_mapping.strip().split(",")])
            source_islands = {
                island.get_id(): island
                for island in read_mge_genomic_islands_gff(args.source_islands)
            }
            dest_islands = {
                island.get_id(): island
                for island in read_mge_genomic_islands_gff(args.dest_islands)
            }
            out_prefix = os.path.join(
                args.output_dir,
                f"{args.genome_id}.mge_islands.liftover"
            )
            gff_out = open(f"{out_prefix}.gff3", "wt", encoding="UTF-8",)

            with gff_out:
                print("##gff-version 3", file=gff_out)
                for id1, dst in dest_islands.items():
                    id2 = island_mapping.get(id1)
                    src = source_islands.get(id2)
                    new_island = GenomicIsland.from_island(dst, dst.genome,)
                    if src is None:
                        new_island.update_recombinases()
                        annotated_island = AnnotatedGenomicIsland.from_island(new_island)
                        mge_island = MgeGenomicIsland.from_island(annotated_island)
                        mge_island.evaluate_recombinases(mge_rules)
                    else:
                        # new_island = GenomicIsland.from_island(dst, dst.genome,)
                        # for src_gene, dst_gene in zip(src.genes, new_island.genes):
                        #     dst_gene.liftover(src_gene)
                        GeneSet.liftover(src.genes, new_island.genes)
                        new_island.update_recombinases()
                        annotated_island = AnnotatedGenomicIsland.from_island(new_island)
                        mge_island = MgeGenomicIsland.from_island(annotated_island)
                        mge_island.evaluate_recombinases(mge_rules)
                    mge_island.to_gff(
                        gff_out,
                        source_db=None,
                    )
                    mge_islands.setdefault(mge_island.contig, []).append(mge_island)
            if args.extract_islands:
                extract_mge_seqs(args.extract_islands, mge_islands, out_prefix)


                


            

        elif args.command == "annotate_genes":
            raise NotImplementedError

        elif args.command == "annotate":
            raise NotImplementedError


if __name__ == "__main__":
    main()
