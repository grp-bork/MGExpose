# pylint: disable=R0913,R0917
""" Module to manage gene sets. """

import logging

from .gene import Gene
from .annotation.eggnog import EMAPPER_FIELDS
from ..utils.gffio import read_prodigal_gff
from ..utils.readers import read_preannotated_genes


logger = logging.getLogger(__name__)


class GeneSet(dict):
    """ Class to manage gene sets. """
    def __init__(
        self,
        genes,
        genome_id=None,
        speci=None,
        composite_gene_ids=False,
    ):
        super().__init__()

        logger.info(
            "Creating new %s for genome=%s specI=%s",
            self.__class__, genome_id, speci,
        )

        for gene in genes:
            if composite_gene_ids:
                # PG3 input is preprocessed (no gffs), so the gene ids are
                # already in the correct format
                # for all other prodigal-based input
                # the gene ids are combined from the contig id and the
                # suffix of col9's ID record:
                # CALOLV020000065.1	[...]	ID=65_14;... -> CALOLV020000065.1_14
                # gene_id = f'{annotation[0]}_{gene_id.split("_")[-1]}'
                # gene.id = f'{gene.contig}_{gene.id.split("_")[-1]}'
                gene.set_composite_id()

            if genome_id is not None and gene.genome is None:
                gene.genome = genome_id
            if speci is not None and gene.speci is None:
                gene.speci = speci

            self[gene.id] = gene

    def dump(self, outstream):
        """ Write gene info to stream. """

        headers = list(Gene().__dict__.keys())
        headers.remove("eggnog")
        headers.remove("secretion_systems")
        headers.remove("secretion_rules")
        headers += ("secretion_systems", "secretion_rules",)
        headers += EMAPPER_FIELDS["v2.1.2"]
        headers.remove("description")

        print(*headers, sep="\t", file=outstream)
        for gene in self.values():
            eggnog_data = {}
            if gene.eggnog:
                eggnog_data = dict(gene.eggnog)
            eggnog_cols = (
                eggnog_data.get(k)
                for k in EMAPPER_FIELDS["v2.1.2"]
                if k != "description"
            )

            secretion_systems, secretion_rules = None, None
            if gene.secretion_systems:
                secretion_systems = ",".join(gene.secretion_systems)
            if gene.secretion_rules:
                secretion_rules = ",".join(str(s) for s in gene.secretion_rules)

            print(
                gene,
                secretion_systems,
                secretion_rules,
                *eggnog_cols,
                sep="\t",
                file=outstream
            )

    @classmethod
    def from_file(
        cls,
        fn,
        genome_id=None,
        speci=None,
        gene_type="prodigal",
        composite_gene_ids=False,
    ):
        """ Generate a GeneSet from a prodigal gff or mgexpose gene_info. """
        if gene_type == "prodigal":
            read_f = read_prodigal_gff
        else:
            read_f = read_preannotated_genes

        return cls(
            read_f(fn),
            genome_id=genome_id,
            speci=speci,
            composite_gene_ids=composite_gene_ids,
        )

    @classmethod
    def from_island(cls, island, composite_gene_ids=False,):
        """ Generate a GeneSet from a GenomicIsland or subclasses. """
        return cls(island.genes, island.genome, island.speci, composite_gene_ids=composite_gene_ids,)
