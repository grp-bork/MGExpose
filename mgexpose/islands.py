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


@dataclass
class AnnotatedGenomicIsland(GenomicIsland):
    '''The following class extends generic Genomic Island with MGE machinery annotations.'''

    TABLE_HEADERS = (
        "contig",
        "start",
        "end",
        "island_size",
        "prot_count",
        "mgeR_count",
        "Plasmid_PA",
        "phage_count",
        "all_conj_count",
        "CONJ_T4SS",
        "SS_present_mandatoryG",
        "entire_ss",
        "mgeR",
    )

    phage_count: int = 0
    conj_count: int = 0
    conj_man_count: int = 0
    valid_entire: bool = False
    valid_mandatory: bool = False
    valid_accessory: bool = False

    def __post_init__(self):
        """ Apply annotations. """
        secretion_systems = {}
        cm_counts = Counter()

        for gene in self.genes:
            self.phage_count += gene.phage is not None
            if gene.secretion_systems:
                self.conj_count += 1

                has_mandatory_system = False
                for system, rule in zip(gene.secretion_systems, gene.secretion_rules):
                    try:
                        _, system = system.split(":")
                    except ValueError:
                        continue
                    if system.split("/")[1].split("_")[0] in ("dCONJ", "T4SS", "MOB",):
                        has_mandatory_system = True
                    if rule is not None:
                        cm_counts[(system, False)] += 1
                        cm_counts[(system, True)] += 1
                        secretion_systems[system] = rule

                self.conj_man_count += has_mandatory_system
                
                logger.info("Secretion system: Gene %s -> conj_count = %s", gene.id, self.conj_man_count)
                # self.conj_man_count += (
                #     # gene.secretion_system.upper()[:4] in ("CONJ", "T4SS",)
                #     gene.secretion_system.split("_")[0] in ("dCONJ", "T4SS", "MOB",)
                # )
                # if gene.secretion_rule is not None:
                #     cm_counts[(gene.secretion_system, False)] += 1
                #     cm_counts[(gene.secretion_system, True)] += 1

                #     secretion_systems[gene.secretion_system] = gene.secretion_rule

        # TODO: validate if still works
        for system, rule in secretion_systems.items():
            self.valid_mandatory |= (cm_counts[(system, True)] >= rule["mandatory"] / 2)
            self.valid_accessory |= (cm_counts[(system, False)] >= rule["accessory"] / 2)
            self.valid_entire |= (
                cm_counts[(system, True)] == rule["mandatory"] and
                cm_counts[(system, False)] == rule["accessory"]
            )

    def __str__(self):
        """ String representation. """
        return "\t".join(
            map(
                str, (
                    self.contig,
                    self.start,
                    self.end,
                    len(self),
                    len(self.genes),
                    # len(self.recombinases),
                    sum(self.recombinases.values()),
                    0,  # is still necessary?
                    self.phage_count,
                    self.conj_count,
                    self.conj_man_count,
                    self.valid_mandatory,
                    self.valid_entire,
                    ",".join(self.recombinases),
                )
            )
        )

    @classmethod
    def from_genomic_island(cls, g_island):
        """ Construct annotated island from genomic island. """
        return cls(
            **g_island.__dict__,
        )


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
        self.recombinases = Counter(recombinases.split(","))

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
        except:
            raise ValueError(f"not enough cols? {cols}")

        try:
            recombinases = Counter(
                dict((key, int(value)) for key, value in
                     (item.split(":")
                      for item in attribs["mgeR"].split(","))
                     )
            )
        except:
            raise ValueError(f"recombinase string weird? {attribs['mgeR'].split(',')}")

        try:
            mges = Counter(
                dict((key, int(value)) for key, value in
                     (item.split(":")
                      for item in attribs["mge"].split(","))
                     )
            )
        except:
            raise ValueError(f"mge string weird? {attribs['mge'].split(',')}")

        if mges.get("is_tn"):
            mges["c_tn"] = mges["is_tn"]
            del mges["is_tn"]

        genome_id, *_ = GenomicIsland.parse_id(attribs["ID"])
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
