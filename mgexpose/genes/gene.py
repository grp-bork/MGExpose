# pylint: disable=R0902,R0917,R0913

""" Gene module """

import re

from ast import literal_eval
from dataclasses import dataclass, field

from .annotation.eggnog import EMAPPER_FIELDS


@dataclass
class Gene:
    '''The following class describes a Gene sequence with its attributes.
    Each gene can contribute to the definition of a MGE Island by being
    1. MGE machinery i.e. phage, conjugation system, conjugation rule
    2. Recombinase i.e. mge (naming is confusing but is kept for historical reasons)
    Eventually each gene has additional annotations coming from EggNOG mapper and
    associated with it and can be extended. '''
    id: str = None
    genome: str = None
    species: str = None
    contig: str = None
    start: int = None
    end: int = None
    strand: str = None
    recombinase: str = None
    cluster: str = None
    is_core: bool = None

    phage: str = None
    eggnog: tuple = None
    conjugation_systems: list = field(default_factory=list)
    conjugation_rules: list = field(default_factory=list)

    parent: str = None

    # specify optional annotations here
    # when adding new class variables,
    # otherwise output will be suppressed.
    OPTIONAL_ANNOTATIONS = (
        "phage",
        "conjugation_systems",
        "conjugation_rules",
        "recombinase",
        "eggnog",
        "parent",
    )
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
        return self.phage is None and self.recombinase is None and not self.conjugation_systems

    def __str__(self):
        """ String representation. """

        species = self.species
        if not isinstance(species, str):
            # converts non-string speci annotation (coreg mode) to string.
            species = ":".join(sorted(species))

        return "\t".join(
            f"{v if k != 'species' else species}" for k, v in self.__dict__.items()
            if k not in ("eggnog", "conjugation_systems", "conjugation_rules",)
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

        conjugation_rules = attribs.get("conjugation_rules")

        genome_type = attribs.get("genome_type")

        return cls(
            id=attribs["ID"],
            genome=attribs.get("genome"),
            species=attribs.get("species", "unknown_species"),
            contig=cols[0],
            start=int(cols[3]),
            end=int(cols[4]),
            strand=cols[6],
            recombinase=attribs.get("recombinase"),
            cluster=attribs.get("cluster") or attribs.get("Cluster"),
            #is_core=(attribs.get("genome_type") == "COR" if attribs.get("genome_type") else None),
            is_core=(genome_type == "COR" if genome_type is not None else genome_type),
            phage=attribs.get("phage"),
            conjugation_systems=attribs.get("conjugation_systems", "").split(","),
            conjugation_rules=literal_eval(f"[{conjugation_rules}]") if conjugation_rules else [],
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
            "conjugation_systems": ",".join(
                self.conjugation_systems
            ) if self.conjugation_systems else None,
            "conjugation_rules": ",".join(
                str(s) for s in self.conjugation_rules
            ) if self.conjugation_rules else None,
            "phage": self.phage,
            "recombinase": self.recombinase,
            "genome_type": self.rtype(self.is_core),
        }
        if intermediate_dump:
            attribs["genome"] = self.genome
            attribs["species"] = self.species
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
        - species
        - contig
        - start
        - end
        - strand
        - recombinase
        - cluster
        - is_core
        - phage
        - conjugation_system
        - conjugation_rule
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

        conjugation_rules = kwargs.get("conjugation_rule", kwargs.get("conjugation_rules"))
        if conjugation_rules is None:
            conjugation_rules = []
        else:
            conjugation_rules = conjugation_rules.strip()
            if re.match(r"\[(\{'mandatory': [0-9]+, 'accessory': [0-9]+\})+\]", conjugation_rules):
                conjugation_rules = literal_eval(conjugation_rules)
            else:
                conjugation_rules = []

        def parse_is_core(s: str):
            if s is None:
                return None
            if not isinstance(s, str):
                raise TypeError(f"{s=} is {type(s)} but has to be string.")
            s = s.capitalize()
            if s not in ("False", "True", "None"):
                return None
            return literal_eval(s)

        # conjugation_systems=attribs.get("conjugation_systems", "").split(","),
        # conjugation_rules=literal_eval(f"[{conjugation_rules}]") if conjugation_rules else [],
        conjugation_systems = kwargs.get(
            "conjugation_system",
            kwargs.get("conjugation_systems", "")
        )
        if conjugation_systems is None:
            conjugation_systems = []

        return cls(
            id=gene_id,
            genome=genome_id,
            species=kwargs.get("species"),
            contig=kwargs.get("contig"),
            start=int(kwargs.get("start")),
            end=int(kwargs.get("end")),
            strand=kwargs.get("strand"),
            recombinase=kwargs.get("recombinase"),
            cluster=kwargs.get("cluster"),
            # is_core=kwargs.get("is_core") == "True",
            is_core=parse_is_core(kwargs.get("is_core", "None")),
            phage=kwargs.get("phage"),
            conjugation_systems=conjugation_systems.split(",") if conjugation_systems else [],
            conjugation_rules=conjugation_rules,
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
