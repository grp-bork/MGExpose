# pylint: disable=R0902,R0917,R0913

""" Gene module """

import copy
import re

from ast import literal_eval
from dataclasses import dataclass, field

from .annotation.eggnog import EMAPPER_FIELDS


@dataclass
class Gene:
    '''The following class describes a Gene sequence with its attributes.
    Each gene can contribute to the definition of a MGE Island by being
    1. MGE machinery i.e. phage, secretion system, secretion rule
    2. Recombinase i.e. mge (naming is confusing but is kept for historical reasons)
    Eventually each gene has additional annotations coming from EggNOG mapper and
    associated with it and can be extended. '''
    id: str = None
    genome: str = None
    speci: str = None
    contig: str = None
    start: int = None
    end: int = None
    strand: str = None
    recombinase: str = None
    cluster: str = None
    is_core: bool = None

    phage: str = None
    eggnog: tuple = None
    secretion_systems: list = field(default_factory=list)
    secretion_rules: list = field(default_factory=list)

    parent: str = None

    # specify optional annotations here
    # when adding new class variables,
    # otherwise output will be suppressed.
    OPTIONAL_ANNOTATIONS = (
        "phage",
        "secretion_systems",
        "secretion_rules",
        "recombinase",
        "eggnog",
        "parent",
    )
    # these are only optional when core genome calculations
    # are disabled, e.g. co-transferred region inputs
    CLUSTER_ANNOTATIONS = ("cluster", "is_core",)

    def liftover(self, other):
        self.recombinase = other.recombinase
        self.cluster = other.cluster
        self.is_core = other.is_core
        self.phage = other.phage
        self.eggnog = copy.deepcopy(other.eggnog) if other.eggnog is not None else None
        self.secretion_systems = copy.deepcopy(other.secretion_systems) if other.secretion_systems is not None else None
        self.secretion_rules = copy.deepcopy(other.secretion_rules) if other.secretion_rules is not None else None

    @staticmethod
    def rtype(is_core):
        """ Returns is_core-tag. """
        if is_core is None:
            return "NA"
        return ("ACC", "COR")[is_core]

    @staticmethod
    def is_core_gene(occ, n_genomes, core_threshold=0.95):
        """ Calculate if a gene is a core gene
        according to its prevalence in a set of species cluster genomes. """
        return occ / n_genomes > core_threshold

    def stringify_eggnog(self):
        """ convert eggnog annotation into gff-col9 key-value pairs """
        if self.eggnog:
            return ";".join(f"{key}={value}" for (key, value) in self.eggnog)
        return None

    def __len__(self):
        """ Calculates gene length. """
        if self.start is None or self.end is None:
            return 0
        return abs(self.end - self.start) + 1

    def is_cargo(self):
        """ Checks if gene can be classified as cargo. """
        return self.phage is None and self.recombinase is None and not self.secretion_systems

    def __str__(self):
        """ String representation. """

        speci = self.speci
        if not isinstance(speci, str):
            # converts non-string speci annotation (coreg mode) to string.
            speci = ":".join(sorted(speci))

        return "\t".join(
            f"{v if k != 'speci' else speci}" for k, v in self.__dict__.items()
            if k not in ("eggnog", "secretion_systems", "secretion_rules",)
        )

    def __hash__(self):
        """ hash function """
        return hash(str(self))

    def has_basic_annotation(self, skip_core_gene_computation=False):
        """ Checks if gene has minimal annotations. """
        ignore = tuple(Gene.OPTIONAL_ANNOTATIONS)
        if skip_core_gene_computation:
            ignore += Gene.CLUSTER_ANNOTATIONS
        for k, v in self.__dict__.items():
            if v is None and k not in ignore:
                return False
        return True

    def is_in_interval(self, start, end):
        """ Checks if gene is located within a region. """
        return start <= self.start <= self.end <= end

    @classmethod
    def from_gff(cls, *cols,):
        """ construct gene from gff record """
        attribs = dict(item.split("=") for item in cols[-1].strip(";").split(";"))

        secretion_rules = attribs.get("secretion_rules")
        return cls(
            id=attribs["ID"],
            genome=attribs.get("genome"),
            speci=attribs.get("speci", "no_speci"),
            contig=cols[0],
            start=int(cols[3]),
            end=int(cols[4]),
            strand=cols[6],
            recombinase=attribs.get("recombinase"),
            cluster=attribs.get("cluster") or attribs.get("Cluster"),
            is_core=attribs.get("genome_type") == "COR",
            phage=attribs.get("phage"),
            secretion_systems=attribs.get("secretion_systems", "").split(","),
            secretion_rules=literal_eval(f"[{secretion_rules}]") if secretion_rules else [],
            eggnog=tuple(
                (k, attribs.get(k))
                for k in EMAPPER_FIELDS["v2.1.2"]
                if attribs.get(k) and k != "description"
            ),
            parent=attribs.get("Parent", "NA",),
        )

    def to_gff(
        self,
        gff_outstream,
        intermediate_dump=False,
        add_header=False,
    ):
        """ dump gene to gff record """

        if add_header:
            print("##gff-version 3", file=gff_outstream)

        attribs = {
            "ID": self.id,
            "Parent": self.parent,
            "cluster": self.cluster,
            "size": len(self),
            "secretion_systems": ",".join(
                self.secretion_systems
            ).strip().strip(",") if self.secretion_systems else None,
            "secretion_rules": ",".join(
                str(s) for s in self.secretion_rules
            ) if self.secretion_rules else None,
            "phage": self.phage,
            "recombinase": self.recombinase,
            "genome_type": self.rtype(self.is_core),
        }
        if intermediate_dump:
            attribs["genome"] = self.genome
            attribs["speci"] = self.speci
            attribs["cluster"] = self.cluster
            attribs["is_core"] = self.is_core

        attrib_str = ";".join(f"{item[0]}={item[1]}" for item in attribs.items() if item[1])

        if self.eggnog:
            attrib_str += f";{self.stringify_eggnog()}"

        print(
            self.contig,
            ".",
            "gene",
            self.start,
            self.end,
            len(self),
            self.strand or ".",
            ".",  # Phase
            attrib_str,
            sep="\t",
            file=gff_outstream,
        )

    @classmethod
    def from_geneinfo(cls, composite_gene_id=False, **kwargs):
        """ Parse a gene from a gene_info line:
        - id
        - genome
        - speci
        - contig
        - start
        - end
        - strand
        - recombinase
        - cluster
        - is_core
        - phage
        - secretion_system
        - secretion_rule
        - cog_fcat
        - seed_eggNOG_ortholog
        - seed_ortholog_evalue
        - seed_ortholog_score
        - eggnog_ogs
        - max_annot_lvl
        - gos
        - ec
        - kegg_ko
        - kegg_pathway
        - kegg_module
        - kegg_reaction
        - kegg_rclass
        - brite
        - cazy
        - bigg_reaction
        - pfam
        """
        gene_id = kwargs.get("id")
        genome_id = kwargs.get("genome")
        if composite_gene_id:
            gene_id = f"{genome_id}.{gene_id}"

        secretion_rules = kwargs.get("secretion_rule", kwargs.get("secretion_rules"))
        if secretion_rules is None:
            secretion_rules = []
        else:
            secretion_rules = secretion_rules.strip()
            if re.match(r"\[(\{'mandatory': [0-9]+, 'accessory': [0-9]+\})+\]", secretion_rules):
                secretion_rules = literal_eval(secretion_rules)
            else:
                secretion_rules = []

        def parse_is_core(s: str):
            if s is None:
                return None
            if not isinstance(s, str):
                raise TypeError(f"{s=} is {type(s)} but has to be string.")
            s = s.capitalize()
            if s not in ("False", "True", "None"):
                return None
            return literal_eval(s)

        # secretion_systems=attribs.get("secretion_systems", "").split(","),
        # secretion_rules=literal_eval(f"[{secretion_rules}]") if secretion_rules else [],
        secretion_systems = kwargs.get("secretion_system", kwargs.get("secretion_systems", ""))
        if secretion_systems is None:
            secretion_systems = []

        return cls(
            id=gene_id,
            genome=genome_id,
            speci=kwargs.get("speci"),
            contig=kwargs.get("contig"),
            start=int(kwargs.get("start")),
            end=int(kwargs.get("end")),
            strand=kwargs.get("strand"),
            recombinase=kwargs.get("recombinase"),
            cluster=kwargs.get("cluster"),
            # is_core=kwargs.get("is_core") == "True",
            is_core=parse_is_core(kwargs.get("is_core", "None")),
            phage=kwargs.get("phage"),
            secretion_systems=secretion_systems.split(",") if secretion_systems else [],
            secretion_rules=secretion_rules,
            eggnog=tuple(
                (k, kwargs.get(k))
                for k in EMAPPER_FIELDS["v2.1.2"]
                if kwargs.get(k) and k != "description"
            ),
            parent=kwargs.get("parent"),
        )

    def set_composite_id(self):
        """
        # PG3 input is preprocessed (no gffs), so the gene ids are
        # already in the correct format
        # for all other prodigal-based input
        # the gene ids are combined from the contig id and the
        # suffix of col9's ID record:
        # CALOLV020000065.1	[...]	ID=65_14;... -> CALOLV020000065.1_14
        """
        # gene_id = f'{annotation[0]}_{gene_id.split("_")[-1]}'
        self.id = f'{self.contig}_{self.id.rsplit("_", maxsplit=1,)[-1]}'
