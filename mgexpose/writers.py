import contextlib
import gzip
import os

from .gffio import write_island_gff
from .islands import MgeGenomicIsland
from .readers import read_fasta


""" module docstring """

MGE_TABLE_HEADERS = \
    ("is_tn",) + \
    MgeGenomicIsland.TABLE_HEADERS[1:6] + \
    MgeGenomicIsland.TABLE_HEADERS[8:14] + \
    ("mgeR", "name", "genes",)



# def dump_islands(islands, out_prefix, db, write_genes=False, add_functional_annotation=False, contig_data=None,):
#     """ dump genomic islands to intermediate gff """
#     with open(
#         f"{out_prefix}.unannotated_islands.gff3",
#         "wt", encoding="UTF-8"
#     ) as _out:
#         write_island_gff(f"{out_prefix}.unannotated_islands.gff3", islands, contig_data, db, write_genes_to_gff=True, add_functional_annotation=True, intermediate_dump=True,)
#         print("##gff-version 3", file=_out)
#         contig = None
#         for island in sorted(islands, key=lambda isl: isl.contig):
#             if island.contig != contig:
#                 contig = island.contig
#                 write_contig(
#                     _out,
#                     contig,
#                     None if contig_data is None else contig_data.get(contig),
#                 )                
#             island.to_gff(
#                 _out, db, write_genes=write_genes,
#                 add_functional_annotation=add_functional_annotation,
#                 intermediate_dump=True,
#             )


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
    contig_data=None,
):
    """ write final results """

    out_prefix = os.path.join(
        output_dir,
        f"{genome_id}.{output_suffix}"
    )

    outstream = contextlib.nullcontext()
    if write_tsv:
        outstream = open(
            f"{out_prefix}.txt",
            "wt",
            encoding="UTF-8",
        )

    islands_by_contig = {}

    with outstream:
        # TSV header
        if write_tsv:
            print(*MGE_TABLE_HEADERS, sep="\t", file=outstream)

        # Start recording the outputs
        # Sort the list of MGEGenomicIslands based on contig names
        for island in sorted(recombinase_islands, key=lambda isl: isl.contig):
            islands_by_contig.setdefault(island.contig, []).append(island)
            # TSV: ignore gene-wise annotations; each line is recombinase island,
            # all gene IDs are stored in a gene_list column
            # assert genome_id == island.genome
            if write_tsv:
                island.to_tsv(outstream)
    
    if write_gff:
        write_island_gff(
            f"{out_prefix}.gff3",
            islands_by_contig,
            contig_data,
            dbformat,
            write_genes_to_gff=write_genes_to_gff,
            add_functional_annotation=add_functional_annotation,
            intermediate_dump=False,
        )
    
    if genome_seqs is not None:
        write_island_sequences(out_prefix, genome_seqs, islands_by_contig)

def write_island_sequences(out_prefix, genome_seqs, islands_by_contig):
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
