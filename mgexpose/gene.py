# pylint: disable=R0902,R0917,R0913

""" Gene module """

from dataclasses import dataclass

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
    secretion_system: str = None
    secretion_rule: dict = None

    # specify optional annotations here
    # when adding new class variables,
    # otherwise output will be suppressed.
    OPTIONAL_ANNOTATIONS = ("phage", "secretion_system", "secretion_rule", "recombinase", "eggnog",)
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
    def is_core_gene(occ, n_genomes, core_threshold=0.95, strict=True):
        if strict or n_genomes == 2 or n_genomes > 20:
            return occ / n_genomes > core_threshold
        return occ >= n_genomes - 1

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
            if k != "eggnog"
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
    def from_gff(cls, *cols):
        """ construct gene from gff record """
        attribs = dict(item.split("=") for item in cols[-1].split(";"))
        return cls(
            id=attribs["ID"],
            genome=attribs.get("genome"),
            speci=attribs.get("speci"),
            contig=cols[0],  # contig
            start=int(cols[3]),  # start
            end=int(cols[4]),  # end
            strand=cols[6],  # strand
            recombinase=attribs.get("recombinase"),
            cluster=attribs.get("cluster") or attribs.get("Cluster"),
            is_core=attribs.get("genome_type") == "COR",
            phage=attribs.get("phage"),
            secretion_system=attribs.get("secretion_system"),
            secretion_rule=attribs.get("secretion_rule"),
            eggnog=tuple(
                (k, attribs.get(k))
                for k in EggnogReader.EMAPPER_FIELDS["v2.1.2"]
                if attribs.get(k) and k != "description"
            ),
        )

    def to_gff(
        self,
        gff_outstream,
        genomic_island_id,
        add_functional_annotation=False,
        intermediate_dump=False,
        add_header=False,
    ):
        """ dump gene to gff record """

        if add_header:
            print("##gff-version 3", file=gff_outstream)

        attribs = {
            "ID": self.id,
            "Parent": genomic_island_id,
            "cluster": self.cluster,
            "size": len(self),
            "secretion_system": self.secretion_system,
            "secretion_rule": self.secretion_rule,
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
