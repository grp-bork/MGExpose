# pylint: disable=R0913,R0914,R0917
""" Module containing writer functions. """

import contextlib
import gzip
import os

from ...islands.mge_genomic_island import MgeGenomicIsland
from ..readers import read_fasta


def dump_islands(islands, out_prefix, db,):
    """ dump genomic islands to intermediate gff """
    with open(
        f"{out_prefix}.unannotated_islands.gff3",
        "wt", encoding="UTF-8"
    ) as _out:
        print("##gff-version 3", file=_out)
        for island in sorted(islands, key=lambda isl: isl.contig):
            island.to_gff(
                _out, db,
                intermediate_dump=True,
            )


def write_final_results(recombinase_islands, args):
    """ write final results """

    outstream = contextlib.nullcontext()

    out_prefix = os.path.join(
        args.output_dir,
        f"{args.genome_id}.mge_islands"
    )

    gff_outstream = open(
        f"{out_prefix}.gff3",
        "wt",
        encoding="UTF-8",
    )

    # Sort the list of MGEGenomicIslands based on contig names
    sorted_islands = sorted(recombinase_islands, key=lambda isl: isl.contig)
    islands_by_contig = {}

    with outstream, gff_outstream:
        # GFF3 header
        print("##gff-version 3", file=gff_outstream)

        # Start recording the outputs
        for island in sorted_islands:
            islands_by_contig.setdefault(island.contig, []).append(island)

            # GFF3: add individual genes annotation;
            # parent lines are recombinase islands, children lines are genes
            # GFF3 parent term: recombinase island
            island.to_gff(
                gff_outstream,
                source_db=args.dbformat,
            )

        if args.genome_fasta:
            extract_mge_seqs(args.genome_fasta, islands_by_contig, out_prefix)


def extract_mge_seqs(genome_seqs, islands, out_prefix):
    """ Extract MGE sequences from genomic sequences. """
    with gzip.open(
        f"{out_prefix}.ffn.gz",
        "wt",
    ) as _out:
        for header, seq in read_fasta(genome_seqs):
            seqid, *_ = header.split(" ")
            for island in islands.get(seqid, []):
                attrib_str = ";".join(
                    f"{item[0]}={item[1]}"
                    for item in island.get_attribs().items()
                    if item[1] and item[0] not in ("ID", "name")
                )
                print(
                    f">{island.get_id()} {attrib_str}",
                    seq[island.start - 1: island.end],
                    sep="\n", file=_out,
                )
