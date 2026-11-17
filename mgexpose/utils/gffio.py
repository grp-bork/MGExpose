""" GFF I/O -- wannabe serialisation module """

from .chunk_reader import get_lines_from_chunks
from ..genes.gene import Gene
from ..islands.genomic_island import GenomicIsland
from ..islands.mge_genomic_island import MgeGenomicIsland


def read_island_gff(fn, island_cls):
    """ Read island gff """
    island = None
    for line in get_lines_from_chunks(fn):
        if line and line[0] != "#":
            cols = line.split("\t")
            if cols[2] in ("region", "mobile_genetic_element"):
                if island is not None:
                    yield island
                island = island_cls.from_gff(*cols)
            elif cols[2] == "gene":
                gene = Gene.from_gff(*cols)
                if island is not None:
                    if not island.genes:
                        island.genes = set()
                    island.add_gene(gene, update_coords=False,)
                    if island.species is None:
                        island.species = gene.species
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
            yield Gene.from_gff(*cols,)
