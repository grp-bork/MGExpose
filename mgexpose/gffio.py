""" GFF I/O -- wannabe serialisation module """

from .gene import Gene
from .islands import GenomicIsland, MgeGenomicIsland
from .utils.chunk_reader import get_lines_from_chunks


def read_island_gff(fn, island_cls):
    """ Read island gff """
    with open(fn, "rt", encoding="UTF-8") as _in:
        island = None
        for line in _in:
            line = line.strip()
            if line and line[0] != "#":
                cols = line.split("\t")
                if cols[2] == island_cls.GFFTYPE:
                    if island is not None:
                        yield island
                    island = island_cls.from_gff(*cols)
                elif cols[2] == "gene":
                    gene = Gene.from_gff(*cols)
                    if island is not None:
                        island.genes.add(gene)
                    else:
                        raise ValueError("Found gene but no island.")
        if island is not None:
            yield island


def read_genomic_islands_gff(fn):
    """ reads a set of genomic islands + genes from a gff3 """
    yield from read_island_gff(fn, GenomicIsland)


def read_mge_genomic_islands_gff(fn):
    """ reads a set of mge genomic islands + genes from a gff3 """
    yield from read_island_gff(fn, MgeGenomicIsland)


def read_prodigal_gff(f):
    """ Prodigal gff output reader.

    Returns Gene objects via generator.
    """
    for line in get_lines_from_chunks(f):
        if line and line[0] != "#":
            cols = line.rstrip().split("\t")
            # _id = [
            #     item.split("=")[1]
            #     for item in line[8].split(";")
            #     if item.startswith("ID")
            # ][0]
            # gene_id = f"{line[0]}_{_id.split('_')[1]}"
            # yield gene_id, line
            # yield _id, line

            yield Gene.from_gff(*cols,)
