#!/usr/bin/env python

# pylint: disable=R0912,R0914,R0915,R0913,R0917

""" Mobile genetic element annotation """

import logging
import os
import pathlib
import sys

from .genes.annotation.annotate import annotate_genes
from .genes.gene_calling import run_pyrodigal
from .genes.geneset import GeneSet
from .handle_args import handle_args
from .islands.annotated_genomic_island import AnnotatedGenomicIsland
from .islands.island_processing import annotate_islands, evaluate_islands, prepare_islands
from .islands.mge_genomic_island import MgeGenomicIsland
from .recombinases.recombinase_scan import run_pyhmmer
from .utils.gffio import read_mge_genomic_islands_gff
from .utils.readers import read_precomputed_islands, read_mge_rules
from .utils.writers import extract_mge_seqs, write_final_results


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

            annotated_genes = annotate_genes(genes, args,)

            precomputed_islands = None
            if args.single_island or args.precomputed_islands:
                precomputed_islands = read_precomputed_islands(
                    genome_id=args.genome_id,
                    single_island=args.single_island,
                    island_file=args.precomputed_islands,
                )

            genomic_islands = prepare_islands(
                annotated_genes, precomputed_islands=precomputed_islands,
            )
            annotated_islands = annotate_islands(genomic_islands,)
            mge_islands = list(evaluate_islands(annotated_islands, read_mge_rules(args.mge_rules),))

            if mge_islands:
                write_final_results(mge_islands, args)

        elif args.command == "reannotate":
            genomic_islands = read_mge_genomic_islands_gff(args.input_gff)
            # genomic_islands = read_island_gff(args.input_gff, GenomicIsland)
            mge_rules = read_mge_rules(args.mge_rules)
            out_prefix = os.path.join(
                args.output_dir,
                f"{args.genome_id}"
            )

            with open(
                f"{out_prefix}.mge_islands.reannotated.gff3", "wt", encoding="UTF-8",
            ) as _out:
                print("##gff-version 3", file=_out)
                mge_islands = {}
                for island in genomic_islands:
                    island = island.to_genomic_island()
                    genes = GeneSet.from_island(island)
                    annotated_genes = annotate_genes(genes, args,)
                    annotated_island = AnnotatedGenomicIsland.from_genomic_island(island)
                    mge_island = MgeGenomicIsland.from_annotated_genomic_island(annotated_island)
                    # mge_island = island
                    mge_island.evaluate_recombinases(mge_rules)
                    mge_island.to_gff(
                        _out,
                        source_db=None,
                        write_genes=True,
                        add_functional_annotation=True,
                    )
                    # mge_islands.append(mge_island)
                    mge_islands.setdefault(mge_island.contig, []).append(mge_island)
                if args.genome_fasta:
                    extract_mge_seqs(args.genome_fasta, mge_islands, out_prefix)

        elif args.command == "annotate_genes":
            raise NotImplementedError

        elif args.command == "annotate":
            raise NotImplementedError


if __name__ == "__main__":
    main()
