# pylint: disable=C0116,C0301,R0902,R0916,R0913,R0917
"""
Data Structures Module

This module is designed to simplify the handling of different genomic sequences,
including but not limited to:

- Genomic Island
- Annotated Genomic Island
- MGE Genomic Island
- Gene

The end product of the pipeline is MGE Genomic Island, an Annotated Genomic Island
consisting of Genes.
It can be saved in a tsv or gff3 format together with its attributes and gene annotations.
The MGE type of each MGE Genomic Island is defined by applying MGE Rule.
"""
import itertools as it
import logging
import sys
import warnings

from collections import Counter
from dataclasses import dataclass, field

from .gene import Gene
from .recombinases import MgeRule, MGE_ALIASES


logger = logging.getLogger(__name__)


@dataclass
class GenomicIsland:
    '''The following class describes a generic genomic region
    with one or more identified recombinases (recombinases).
    This region is then referred as Recombinase Island.
    The Genomic Island is represented by contig, start and end
    coordinates, set of genes, some of which are recombinases and MGE machinery.
    Importantly, the set of genes does not include the non-coding regions.
    '''
    RAW_TABLE_HEADER = (
        "specI",
        "genome_accession",
        "panG",
        "contig",
        "start",
        "end",
        "gene_list",
    )
    GFFTYPE = "region"

    speci: str = None
    genome: str = None
    is_core: bool = None
    contig: str = None
    start: int = None
    end: int = None
    name: str = "ISLAND"

    genes: set = field(default_factory=set)
    # recombinases: list = field(default_factory=list)
    recombinases: Counter = field(default_factory=Counter)

    @staticmethod
    def parse_id(id_string):
        """ Parse genome id, contig id, start and end coordinates from id string.
         Reverses get_id(). """
        cols = id_string.split("_")
        contig, coords = cols[3].split(':')

        return "_".join(cols[1:3]), contig, int(coords[0]), int(coords[1])

    @staticmethod
    def get_fieldnames():
        """ Returns column headers for island table. """
        return (
            "first_recombinase_gene",
            "first_recombinase",
            "island_size",
            "genome",
            "specI",
            "core_acc",
            "contig",
            "first_gene_start",
            "last_gene_end",
            "protid_gene_clusters",
        )

    @classmethod
    def from_region_string(cls, region, genome_id=None,):
        """ Creates island from a predefined region string. """
        _, _, contig, start_end, *_ = region.strip().split(".")
        contig = contig.split(".")[-1]
        start, end = map(int, start_end.split("-"))
        return cls(None, genome_id, None, contig, start, end, region.strip())

    @classmethod
    def from_gene(cls, gene):
        """ Creates island from starting gene. """
        island = cls(
            gene.speci,
            gene.genome,
            gene.is_core,
            gene.contig,
            gene.start,
            gene.end,
        )
        island.add_gene(gene)
        return island

    def __len__(self):
        """ Calculates island length. """
        if self.start is None or self.end is None:
            return 0
        return abs(self.end - self.start) + 1

    def __str__(self):
        """ String representation. """
        genes = (
            f"{gene.id}.{gene.cluster}"
            for gene in sorted(
                self.genes, key=lambda g: (g.start, g.end, g.strand)
            )
        )

        return "\t".join(
            [
                f"{v}" if (k != "is_core" or v is None) else Gene.rtype(self.is_core)
                for k, v in self.__dict__.items()
                if k not in ("genes", "recombinases")
            ] + [",".join(genes)]
        )

    def add_gene(self, gene):
        """ Adds a gene to the island. """
        if gene not in self.genes:
            self.end = max(self.end, gene.end)
            if gene.recombinase:
                # self.recombinases.append(
                #     (f"{gene.id}.{gene.cluster}", gene.recombinase)
                # )
                self.recombinases[gene.recombinase] += 1
            self.genes.add(gene)

    def get_position(self):
        """ Return genomic position tuple. """
        return (self.contig, self.start, self.end)

    def get_recombinases(self):
        for g in sorted(self.genes, key=lambda x: x.start):
            if g.recombinase:
                yield f"{g.id}.{g.cluster}", g.recombinase

    def dump(self, seen_islands, raw_outstream=None, outstream=sys.stdout):
        """ Writes island to outstream. """
        if raw_outstream:
            print(self, file=raw_outstream)
            pos = self.get_position()
            if pos not in seen_islands and self.recombinases:
                seen_islands.add(pos)

                print(
                    # *self.recombinases[0],
                    *tuple(self.get_recombinases())[0],
                    len(self),
                    str(self),
                    sep="\t",
                    file=outstream,
                )

    def get_id(self):
        return f"GIL_{self.genome}_{self.contig}:{self.start}-{self.end}"

    @classmethod
    def from_gff(cls, *cols):
        attribs = dict(item.split("=") for item in cols[-1].split(";"))
        recombinases = Counter(
            {
                item.split(":")[0]: int(item.split(":")[1])
                for item in attribs["recombinases"].split(",")
            }
        )

        return cls(
            attribs["specI"],
            attribs["genome"],
            attribs["genome_type"] == "COR",
            cols[0],  # contig
            int(cols[3]),  # start
            int(cols[4]),  # end
            genes=set(),
            recombinases=recombinases,
        )

    def to_gff(
        self,
        gff_outstream,
        source_db,
        write_genes=False,
        add_functional_annotation=False,
        intermediate_dump=False,
        add_header=False,
    ):

        if add_header:
            print("##gff-version 3", file=gff_outstream)

        island_id = self.get_id()
        attribs = {
            "ID": island_id,
            "genome": self.genome,
            "genome_type": Gene.rtype(self.is_core),
            "size": len(self),
            "n_genes": len(self.genes),
            # "mgeR": ",".join(sorted(r for _, r in self.recombinases)),
            # "mgeR": ",".join(sorted(self.recombinases)),
            "recombinases": (
                ",".join(
                    f"{k}:{v}"
                    for k, v in sorted(self.recombinases.items())
                )
                if self.recombinases else ""
            ),
            "specI": self.speci,
        }
        if self.name:
            attribs["name"] = self.name
        attrib_str = ";".join(f"{item[0]}={item[1]}" for item in attribs.items() if item[1])
        # Format the source column
        if source_db:
            source = f"proMGE_{source_db}"
        else:
            source = "proMGE"
        print(
            self.contig,
            source,
            GenomicIsland.GFFTYPE,
            self.start,
            self.end,
            len(self),  # Score field
            ".",  # Strand
            ".",  # Phase
            attrib_str,
            sep="\t",
            file=gff_outstream
        )

        if write_genes:
            # GFF3 child term: genes
            for gene in sorted(self.genes, key=lambda g: (g.start, g.end,)):
                gene.to_gff(
                    gff_outstream,
                    # genomic_island_id=island_id,
                    add_functional_annotation=add_functional_annotation,
                    intermediate_dump=intermediate_dump,
                )
