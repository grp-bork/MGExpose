# pylint: disable=C0116,C0301,R0902,R0916,R0913,R0917,W0613
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
import warnings

from collections import Counter
from dataclasses import dataclass

from .annotated_genomic_island import AnnotatedGenomicIsland
from .genomic_island import GenomicIsland
from ..genes.gene import Gene
from ..recombinases import MgeRule, MGE_ALIASES


logger = logging.getLogger(__name__)


@dataclass
class MgeGenomicIsland(AnnotatedGenomicIsland):
    '''The following class describes Anotated Genomic Islands with their assigned MGE type.
    Those are Mobile Genetic Elements (MGEs).
    The class attributes are used to describe the MGE properties.
    It also contains functionality to save the MGEs in gff3 or tsv formats.'''

    TABLE_HEADERS = (
        "tn",
        "phage",
        "phage_like",
        "ce",
        "integron",
        "mi",
        "nmi",
        "nov",
        "cellular",
        "contig",
        "start",
        "end",
        "size",
        "n_genes",
        "phage_count",
        "conj_man_count",
        "recombinases",
    )
    GFFTYPE = "mobile_genetic_element"

    integron: int = 0
    cellular: int = 0
    phage: int = 0
    c_mi: int = 0
    nov: int = 0
    c_pli: int = 0
    c_ce: int = 0
    c_nmi: int = 0
    c_tn: int = 0

    tn3_found: bool = False
    ser_found: bool = False

    def __post_init__(self):
        """ Apply annotations. """
        recombinases = (",".join(it.chain(*((r,) * v for r, v in self.recombinases.items())))).lower()
        for name, alias in MGE_ALIASES.items():
            recombinases = recombinases.replace(name, alias)

        self.tn3_found = "tn3" in recombinases
        self.ser_found = "c2_n1ser" in recombinases or "ser_ce" in recombinases

        # integron
        self.integron = int("integron" in recombinases)

        # self.recombinases = recombinases.split(",") if recombinases else []
        self.recombinases.clear()
        self.recombinases.update(recombinases.split(","))
        # self.recombinases = Counter(recombinases.split(","))

        # tag recombinase island with more than 3 recombinases
        self.c_nmi = sum(self.recombinases.values()) > 3

    def __str__(self):
        """ String representation. """
        return "\t".join(
            tuple(map(str, self.get_mge_metrics())) +
            (
                self.contig,
                f"{self.start}",
                f"{self.end}",
                f"{len(self)}",
                f"{len(self.genes)}",
                f"{self.phage_count}",
                f"{self.conj_man_count}",
                # ",".join(self.recombinases),
                ",".join(
                    f"{k}:{v}"
                    for k, v in sorted(self.recombinases.items())
                )
                if self.recombinases else "",
                self.name,
            )
        )

    def get_mge_metrics(self):
        """ Cast mge metrics to int. """
        return tuple(
            map(
                int,
                (
                    self.c_tn,
                    self.phage,
                    self.c_pli,
                    self.c_ce,
                    self.integron,
                    self.c_mi,
                    self.cellular,
                )
            )
        )

    def get_annotated_mge_metrics(self):
        metrics = list(self.get_mge_metrics())  # Get mge_type and counts
        mge_metrics = [
            (k, v)
            for k, v in zip(
                ("is_tn", "phage", "phage_like", "ce", "integron", "mi", "cellular",),
                metrics
            )
            if v  # Collect as long as metrics are not None
        ]
        return mge_metrics

    @staticmethod
    def is_nested(annotated_mge_metrics):
        n_mges = sum(v for _, v in annotated_mge_metrics)
        if not n_mges:
            # raise UserWarning("No MGEs were assigned to recombinase island")
            warnings.warn("No MGEs were assigned to recombinase island")
        # Solitary or nested MGE?
        return n_mges > 1

    @staticmethod
    def mge_num_island_type(is_nested):
        """ Returns nested vs solitary MGE-tag. """
        return ("non-nested", "nested")[is_nested]

    def has_annotation(self):
        """ (Sanity) Check if island has any mge annotation. """
        return sum((
            self.c_tn,
            self.phage,
            self.c_pli,
            self.c_ce,
            self.integron,
            self.c_mi,
            self.cellular,
        )) > 0

    def evaluate_recombinases(self, rules, outstream=None, outstream2=None):
        """ Annotate recombinases. """
        patch_c_tn = False

        recombinases = it.chain(*it.chain((r,) * c for r, c in self.recombinases.items()))

        for rec in recombinases:
            rule = rules.get(rec)
            if rule is None:
                print(f"WARNING: cannot find mge-rule for `{rec}`")
                rule = MgeRule()

            # cellular:Arch1/Cyan/Xer/Cand
            self.cellular |= rule.cellular

            self.c_tn = rule.c_tn_check(self)
            patch_c_tn |= rule.patch_c_tn_check(self)

            if self.phage_count >= 2 and self.conj_man_count < 1:
                self.phage, self.c_mi, self.nov = rule.phage_check(self)
            elif self.phage_count < 2 and self.conj_man_count < 1:
                self.c_pli, self.c_mi = rule.phage_like_check(
                    self,
                    "brujita" in rec
                )
            elif self.phage_count < 2 and self.conj_man_count >= 1:
                self.c_ce, self.nov = rule.conjug_element_check(self)
            elif self.phage_count >= 2 and self.conj_man_count >= 1:
                self.phage, self.c_mi, self.nov = rule.mobility_island_check(self)

        # counting multiple tn in Tn3 containing recombinase island
        # self.c_tn += (len(self.recombinases) > 2) * (self.tn3_found or self.ser_found)
        self.c_tn += (sum(self.recombinases.values()) > 2) * (self.tn3_found or self.ser_found)
        if not self.has_annotation():
            if not patch_c_tn:
                print(f"WARNING: No annotation found, but cannot patch either.\n{self}")
            self.c_tn = patch_c_tn

        if outstream:
            print(self, sep="\t", file=outstream,)

        # previous step in some cases generates overlap between Phage/Phage_like and Mobility island
        # this step specifically resolves such instances based on recombinase presence and presence/
        # absence of phage structural genes/conjugation machinery genes in the neighbourhood
        if self.c_mi and self.c_pli:
            self.c_mi = int(
                any(
                    pat in ",".join(self.recombinases).lower()
                    for pat in ('relaxase', 'rep_', 'mob_', 'trwc')
                )
            )
            self.c_pli = int(not self.c_mi)

        if self.phage and self.c_mi and self.phage_count >= 2:
            self.phage, self.c_mi = True, False

        if outstream2:
            print(self, sep="\t", file=outstream2,)

    @classmethod
    def from_annotated_genomic_island(cls, ag_island):
        """ Construct from annotated genomic island. """
        island = cls(
            **ag_island.__dict__
        )
        for gene in island.genes:
            gene.parent = island.get_id()
        return island

    def get_id(self):
        return f"MGE_{self.genome}_{self.contig}:{self.start}-{self.end}"

    def get_attribs(self):
        mge_metrics = self.get_annotated_mge_metrics()
        attribs = {
            "ID": self.get_id(),
            "mge": ",".join(f"{k}:{v}" for k, v in mge_metrics),  # Count each mge type
            "genome_type": Gene.rtype(self.is_core),
            "mge_type": self.mge_num_island_type(self.is_nested(mge_metrics)),
            "size": len(self),
            "n_genes": len(self.genes),
            "mgeR": (
                ",".join(
                    f"{k}:{v}"
                    # for k, v in sorted(Counter(self.recombinases).items())
                    for k, v in sorted(self.recombinases.items())
                )
                if self.recombinases else ""
            ),

        }
        if self.name:
            attribs["name"] = self.name
        return attribs

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

        # island_id = self.get_id()
        # mge_metrics = self.get_annotated_mge_metrics()
        # attribs = {
        #     "ID": island_id,
        #     "mge": ",".join(f"{k}:{v}" for k, v in mge_metrics),  # Count each mge type
        #     "genome_type": Gene.rtype(self.is_core),
        #     "mge_type": self.mge_num_island_type(self.is_nested(mge_metrics)),
        #     "size": len(self),
        #     "n_genes": len(self.genes),
        #     "mgeR": (
        #         ",".join(
        #             f"{k}:{v}"
        #             # for k, v in sorted(Counter(self.recombinases).items())
        #             for k, v in sorted(self.recombinases.items())
        #         )
        #         if self.recombinases else ""
        #     ),
        # }
        # if self.name:
        #     attribs["name"] = self.name
        attribs = self.get_attribs()
        attrib_str = ";".join(f"{item[0]}={item[1]}" for item in attribs.items() if item[1])
        # Format the source column
        source = ("proMGE", f"proMGE_{source_db}")[bool(source_db)]
        # if source_db:
        #     source = f"proMGE_{source_db}"
        # else:
        #     source = "proMGE"
        print(
            self.contig,
            source,
            MgeGenomicIsland.GFFTYPE,
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
                    # genomic_island_id=attribs["ID"],
                    add_functional_annotation=add_functional_annotation,
                )

    @classmethod
    def from_gff(cls, *cols):
        try:
            attribs = dict(item.split("=") for item in cols[-1].split(";"))
        except Exception as exc:
            raise ValueError(f"not enough cols? {cols}") from exc

        try:
            recombinases = Counter(
                dict((key, int(value)) for key, value in
                     (item.split(":")
                      for item in attribs["mgeR"].split(","))
                     )
            )
        except Exception as exc:
            raise ValueError(f"recombinase string weird? {attribs['mgeR'].split(',')}") from exc

        try:
            mges = Counter(
                dict((key, int(value)) for key, value in
                     (item.split(":")
                      for item in attribs["mge"].split(","))
                     )
            )
        except Exception as exc:
            raise ValueError(f"mge string weird? {attribs['mge'].split(',')}") from exc

        # print(mges)
        for mge_type, mge_type_ in [
            ("is_tn", "c_tn"), ("mi", "c_mi"), ("phage_like", "c_pli"),
            ("ce", "c_ce"), ("nmi", "c_nmi"),
        ]:
            if mges.get(mge_type):
                mges[mge_type_] = mges[mge_type]
                del mges[mge_type]

        genome_id, *_ = GenomicIsland.parse_id(attribs["ID"], cols[0])
        # TODO: check coordinates and ID overlap
        return cls(
            "",  # TODO: where to get/ how to handle specI
            genome_id,
            attribs["genome_type"] == "COR",
            cols[0],  # contig
            int(cols[3]),  # start
            int(cols[4]),  # end
            recombinases=recombinases,
            # mge=mges,
            **mges,
            # mge_type=attribs["mge_type"],
            # size=int(attribs["size"]),
            # n_genes=int(attribs["n_genes"]),
            genes=set(),
        )

    def to_tsv(self, outstream):
        metrics = list(self.get_mge_metrics())
        print(
            *metrics,
            self.contig,
            self.start,
            self.end,
            len(self),  # size
            len(self.genes),  # n_genes
            ",".join(
                f"{k}:{v}"
                # for k, v in sorted(Counter(self.recombinases).items())
                for k, v in sorted(self.recombinases.items())
            ) if self.recombinases else "",
            (self.name if self.name else ""),
            ",".join(gene.id for gene in sorted(self.genes, key=lambda g: g.id)),  # gene_list
            sep="\t",
            file=outstream,
        )

    def to_genomic_island(self):
        return GenomicIsland(
            self.speci,
            self.genome,
            self.is_core,
            self.contig,
            self.start,
            self.end,
            self.name,
            self.genes,
            self.recombinases,
        )
