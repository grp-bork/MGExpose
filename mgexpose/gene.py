# pylint: disable=R0902,R0917,R0913

""" Gene module """

import re

from ast import literal_eval
from dataclasses import dataclass, field

from .eggnog import EggnogReader


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
    OPTIONAL_ANNOTATIONS = ("phage", "secretion_systems", "secretion_rules", "recombinase", "eggnog", "parent",)
    # these are only optional when core genome calculations
    # are disabled, e.g. co-transferred region inputs
    CLUSTER_ANNOTATIONS = ("cluster", "is_core",)

    @staticmethod
    def rtype(is_core):
        """ Returns is_core-tag. """
        if is_core is None:
            return "NA"
        return ("ACC", "COR")[is_core]

    @staticmethod
    def is_core_gene(occ, n_genomes, core_threshold=0.95):
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

    def __str__(self):
        """ String representation. """
        return "\t".join(
            f"{v}" for k, v in self.__dict__.items()
            if k not in ("eggnog", "secretion_systems", "secretion_rules",)
        )

    def stringify_speci(self):
        """ Converts non-string speci annotation (coreg mode) to string. """
        if not isinstance(self.speci, str):
            self.speci = ":".join(sorted(self.speci))

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
            speci=attribs.get("speci"),
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
                for k in EggnogReader.EMAPPER_FIELDS["v2.1.2"]
                if attribs.get(k) and k != "description"
            ),
            parent=attribs.get("Parent", "NA",),
        )

    def to_gff(
        self,
        gff_outstream,
        add_functional_annotation=False,
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
            "secretion_systems": ",".join(self.secretion_systems) if self.secretion_systems else None,
            "secretion_rules": ",".join(str(s) for s in self.secretion_rules) if self.secretion_rules else None,
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

        if add_functional_annotation and self.eggnog:
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
        # id      genome  speci   contig  start   end     strand  recombinase     cluster is_core phage   secretion_system        secretion_rule  cog_fcat        seed_eggNOG_ortholog    seed_ortholog_evalue    seed_ortholog_score     eggnog_ogs      max_annot_lvl   goes    ec      kegg_ko kegg_pathway    kegg_module     kegg_reaction   kegg_rclass     brite   cazy    bigg_reaction   pfam
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
            if not isinstance(s, str):
                raise TypeError(f"{s=} is {type(s)} but has to be string.")
            if s is None or s.lower().capitalize() not in ("False", "True", "None"):
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
                for k in EggnogReader.EMAPPER_FIELDS["v2.1.2"]
                if kwargs.get(k) and k != "description"
            ),
            parent=kwargs.get("parent"),
        )
