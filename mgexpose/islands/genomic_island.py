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

import logging
import re
import sys

from collections import Counter
from dataclasses import dataclass, field

from ..genes.gene import Gene
from ..recombinases import parse_recombinase_string
from ..utils.gffutils import get_attrib_str, get_source_column, parse_gff_attribs


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
    ID_PREFIX = "GIL"


    speci: str = None
    genome: str = None
    is_core: bool = None
    contig: str = None
    start: int = None
    end: int = None
    name: str = "ISLAND"

    genes: set = field(default_factory=set)
    recombinases: Counter = field(default_factory=Counter)

    @staticmethod
    def parse_id(id_string, contig_id=None,):
        """ Parse genome id, contig id, start and end coordinates from id string.
         Reverses get_id(). """
        
        island_id_match = re.match(r'(GIL|MGE|SPIRE)_(.+)_(.+):(\d+)-(\d+)', id_string)
        if not island_id_match:
            raise ValueError(f"{id_string} does not seem to be a valid island identifier.")
        
        groups = island_id_match.groups()[1:]

        start, end = map(int, groups[-2:])

        if contig_id:
            genome_id = "_".join(groups[:2]).replace(f"_{contig_id}", "")
        else:
            genome_id, contig_id = groups[:2]

        return genome_id, contig_id, start, end    


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

    def add_gene(self, gene, update_coords=True,):
        """ Adds a gene to the island. """
        if gene not in self.genes:
            if update_coords:
                self.end = max(self.end, gene.end)
            if gene.recombinase:
                self.recombinases[gene.recombinase] += 1
            self.genes.add(gene)

    def update_recombinases(self):
        self.recombinases.update(
            gene.recombinase
            for gene in self.genes
            if gene.recombinase
        )

    def has_compatible_genome_type(self, gene):
        return self.is_core == gene.is_core

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
        return f"{self.ID_PREFIX}_{self.genome}_{self.contig}:{self.start}-{self.end}"
    
    def get_attribs(self):
        attribs = {
            "ID": self.get_id(),
            "genome": self.genome,
            "genome_type": Gene.rtype(self.is_core),
            "size": len(self),
            "n_genes": len(self.genes),
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

        return attribs


    @classmethod
    def from_gff(cls, *cols):
        attribs = parse_gff_attribs(cols[-1])
        recombinases = attribs.get("recombinases", "")
        recombinases = parse_recombinase_string(recombinases) if recombinases else Counter()
        
        return cls(
            attribs.get("specI"),
            attribs.get("genome"),
            attribs.get("genome_type") == "COR",
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
        intermediate_dump=False,
        add_header=False,
    ):

        if add_header:
            print("##gff-version 3", file=gff_outstream)        

        print(
            self.contig,
            get_source_column(source_db=source_db,),
            GenomicIsland.GFFTYPE,
            self.start,
            self.end,
            len(self),  # Score field
            ".",  # Strand
            ".",  # Phase
            get_attrib_str(self.get_attribs()),
            sep="\t",
            file=gff_outstream,
        )

        for gene in sorted(self.genes, key=lambda g: (g.start, g.end,)):
            gene.to_gff(
                gff_outstream,
                intermediate_dump=intermediate_dump,
            )


    @classmethod
    def from_island(cls, other, genome_id=None, set_parent=True,):
        island = cls(
            **{
                k: other.__dict__.get(k)
                for k in set(other.__dataclass_fields__).intersection(cls.__dataclass_fields__)
            }
        )
        if genome_id:
            island.genome = genome_id
        if set_parent:
            for gene in island.genes:
                gene.parent = island.get_id()
                gene.genome = island.genome
        # for gene in island.genes:
        #     if gene.recombinase:
        #         island.recombinases[gene.recombinase] += 1

        return island
