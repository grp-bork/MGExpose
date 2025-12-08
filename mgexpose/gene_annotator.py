# pylint: disable=R0912,R0913,R0914,R0917

""" Classes for integrating gene annotations. """

import logging
import sys

from contextlib import nullcontext

from .clustering_parser import parse_full_seq_clusters, parse_y_clusters, parse_db_clusters, evaluate_y_clusters
from .gene import Gene
from .phage import PhageDetection
from .eggnog import EggnogReader
from .readers import (
    parse_macsyfinder_report,
    read_recombinase_hits,
)


logger = logging.getLogger(__name__)


class GeneAnnotator:
    """ GeneAnnotator class. """
    def __init__(
        self,
        genome_id,
        speci,
        genes,
        include_genome_id=False,
        has_batch_data=False,
        composite_gene_ids=False,
    ):
        logger.info("Creating new %s for genome=%s specI=%s", self.__class__, genome_id, speci)
        self.genome_id = genome_id
        self.speci = speci
        self.genes = {}
        self.has_batch_data = has_batch_data
        self.include_genome_id = include_genome_id

        # for gene_id, annotation in genes:
        for gene in genes:

            if composite_gene_ids:
                # PG3 input is preprocessed (no gffs), so the gene ids are
                # already in the correct format
                # for all other prodigal-based input
                # the gene ids are combined from the contig id and the
                # suffix of col9's ID record:
                # CALOLV020000065.1	[...]	ID=65_14;... -> CALOLV020000065.1_14
                # gene_id = f'{annotation[0]}_{gene_id.split("_")[-1]}'
                gene.id = f'{gene.contig}_{gene.id.split("_")[-1]}'

            logger.info("Adding gene %s", gene.id)

            gene.genome = self.genome_id
            if gene.speci is None:
                gene.speci =  self.speci
            self.genes[gene.id] = gene

            # self.genes[gene_id] = Gene(
            #     id=gene_id,
            #     genome=self.genome_id,
            #     speci=self.speci,
            #     contig=annotation[0],
            #     start=int(annotation[3]),
            #     end=int(annotation[4]),
            #     strand=annotation[6],
            # )

    def add_recombinases(self, recombinases):
        """ Add information from recombinase scan """
        for gene_id, recombinase in recombinases:
            gene = self.genes.get(gene_id)
            if gene is not None:
                gene.recombinase = recombinase

    def add_cluster(
        self,
        cluster_data,
        use_y_clusters=False,
        core_threshold=0.95,
        output_dir=None,
    ):
        """ Add information from gene clustering to allow for core/accessory gene classification """

        if use_y_clusters:
            if core_threshold == -1:
                parse_y_clusters(cluster_data, self.genes)
            else:
                evaluate_y_clusters(cluster_data, self.genes, core_threshold=core_threshold,)
            return None

        write_data = False
        gene_clusters_out = nullcontext()
        n_genomes = 0
        cluster_genes = {}

        with gene_clusters_out:
            if cluster_data is not None:

                if core_threshold != -1:
                    n_genomes, cluster_genes, gene_cluster_map, _ = parse_full_seq_clusters(
                        self.genome_id,
                        self.genes,
                        cluster_data,
                        output_dir=output_dir,
                    )

                    logger.info(
                        "Parsed %s genomes with %s gene clusters.",
                        n_genomes,
                        len(cluster_genes),
                    )
                else:
                    gene_cluster_map = parse_db_clusters(cluster_data)

                    n_genes = len(gene_cluster_map)
                    n_core_genes = sum(1 for _, _, is_core in gene_cluster_map if is_core)
                    logger.info(
                        "Parsed %s precomputed gene-cluster mappings with %s core genes (%s%%)",
                        n_genes,
                        n_core_genes,
                        round(n_core_genes / n_genes, 2),
                    )

                for gene_id, *cluster in gene_cluster_map:
                    cluster, *is_core = cluster
                    is_core = is_core[0].lower() == "true" if is_core else None
                    if not self.include_genome_id or gene_id.startswith(self.genome_id):
                        gene = self.genes.get(
                            gene_id,
                            self.genes.get(
                                gene_id.replace(self.genome_id + ".", "")
                            )
                        )
                        logger.info(
                            "Checking cluster %s gene %s... %s",
                            str(cluster),
                            gene_id,
                            str(gene),
                        )
                        if gene and gene.speci is not None:
                            gene.cluster = cluster

                            if cluster_genes:
                                occ = cluster_genes[cluster]
                                gene.is_core = Gene.is_core_gene(occ, n_genomes, core_threshold=core_threshold,)
                            elif core_threshold == -1:
                                gene.is_core = is_core

                            if write_data:
                                print(gene.id, gene.cluster, sep="\t", file=gene_clusters_out)

        return None

    def add_eggnog_annotation(self, eggnog_annotation):
        """ Add eggnog output and phage signals to each gene """
        for gene_id, phage_data, eggnog_data in eggnog_annotation:
            gene = self.genes.get(gene_id)
            if gene is not None:
                gene.eggnog = eggnog_data
                gene.phage = phage_data

    def add_secretion_systems(self, secretion_annotation):
        """ Add information from txsscan """
        for gene_id, secretion_data in secretion_annotation:
            
            gene = self.genes.get(gene_id)
            if gene is not None:
                for sgene, system, rule, *_ in secretion_data:
                    gene.secretion_systems.append(f"{sgene}:{system}")
                    gene.secretion_rules.append(rule)

    def annotate_genes(
            self,
            recombinases,
            eggnog_annotation,
            secretion_annotation,
            clusters=None,
            use_y_clusters=False,
            core_threshold=None,
            output_dir=None,
            pyhmmer=True,
    ):
        """ Annotate genes with MGE-relevant data. """
        if recombinases is not None:
            self.add_recombinases(
                read_recombinase_hits(recombinases, pyhmmer=pyhmmer,)
            )
        if secretion_annotation is not None and all(secretion_annotation):
            self.add_secretion_systems(
                parse_macsyfinder_report(
                    *secretion_annotation[:2],
                    # macsy_version=secretion_annotation[-1],
                ),
            )
        if eggnog_annotation is not None:
            phage_detection = PhageDetection(eggnog_annotation[1])

            self.add_eggnog_annotation(
                EggnogReader.parse_emapper(
                    eggnog_annotation[0],
                    phage_annotation=phage_detection,
                )
            )
        if clusters is not None:
            self.add_cluster(
                clusters,
                use_y_clusters=use_y_clusters,
                core_threshold=core_threshold,
                output_dir=output_dir,
            )
        # yield from self.genes.values()
        return self.genes.values()

    def dump_genes(self, outstream):
        """ Write gene info to stream. """

        headers = list(Gene().__dict__.keys())
        headers.remove("eggnog")
        headers.remove("secretion_systems")
        headers.remove("secretion_rules")
        headers += ("secretion_systems", "secretion_rules",)
        headers += EggnogReader.EMAPPER_FIELDS["v2.1.2"]
        headers.remove("description")

        print(*headers, sep="\t", file=outstream)
        for gene in self.genes.values():
            gene.stringify_speci()
            eggnog_data = {}
            if gene.eggnog:
                eggnog_data = dict(gene.eggnog)
            eggnog_cols = (
                eggnog_data.get(k)
                for k in EggnogReader.EMAPPER_FIELDS["v2.1.2"]
                if k != "description"
            )

            secretion_systems = ",".join(gene.secretion_systems) if gene.secretion_systems else None
            secretion_rules = ",".join(str(s) for s in gene.secretion_rules) if gene.secretion_rules else None

            print(gene, secretion_systems, secretion_rules, *eggnog_cols, sep="\t", file=outstream)
