# pylint: disable=R0913,R0914,R0917
""" Module containing writer functions. """

import contextlib
import gzip
import os

from ...islands.mge_genomic_island import MgeGenomicIsland
from ..readers import read_fasta


MGE_TABLE_HEADERS = \
    ("is_tn",) + \
    MgeGenomicIsland.TABLE_HEADERS[1:6] + \
    MgeGenomicIsland.TABLE_HEADERS[8:14] + \
    ("mgeR", "name", "genes",)


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
        f"{args.genome_id}.{args.output_suffix}"
    )

    write_tsv = True
    if write_tsv:
        outstream = open(
            f"{out_prefix}.txt",
            "wt",
            encoding="UTF-8",
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
        # TSV header
        if write_tsv:
            print(*MGE_TABLE_HEADERS, sep="\t", file=outstream)
        # GFF3 header
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
            island.to_gff(
                gff_outstream,
                source_db=args.dbformat,
            )

        genome_seqs = args.extract_islands
        if genome_seqs is not None:
            extract_mge_seqs(genome_seqs, islands_by_contig, args.out_prefix)


def extract_mge_seqs(genome_seqs, islands, out_prefix):
    """ Extract MGE sequences from genomic sequences. """
    with gzip.open(
        f"{out_prefix}.ffn.gz",
        "wt",
    ) as _out:
        for header, seq in read_fasta(genome_seqs):
            seqid, *_ = header.split(" ")
            for island in islands.get(seqid, []):
                attribs = island.get_attribs()
                try:
                    del attribs["ID"]
                except KeyError:
                    pass
                try:
                    del attribs["name"]
                except KeyError:
                    pass
                attrib_str = ";".join(
                    f"{item[0]}={item[1]}" for item in attribs.items() if item[1]
                )
                print(
                    f">{island.get_id()} {attrib_str}",
                    seq[island.start - 1: island.end],
                    sep="\n", file=_out,
                )
