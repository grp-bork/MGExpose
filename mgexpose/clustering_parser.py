# pylint: disable=R0902,R0914

""" Functions for gene cluster parsing """

import logging
import os

from collections import Counter
from contextlib import nullcontext
from dataclasses import dataclass

from .utils.chunk_reader import get_lines_from_chunks


logger = logging.getLogger(__name__)

def evaluate_y_clusters(cluster_data, genes, core_threshold=0.95,):
    # print("EVALUATE_Y_CLUSTERS")
    # print(*list(genes.items())[:10], sep="\n")
    for line in get_lines_from_chunks(cluster_data):
        last_cluster, float_frac = None, None
        gene_id, cluster, _, _, _, frac = line.strip().split("\t")
        if cluster != last_cluster:
            last_cluster, float_frac = cluster, float(frac)
        gene = genes.get(gene_id)
        # print(gene_id, cluster, frac, "->", gene)
        if gene is not None:
            gene.cluster = f"{cluster}:{frac}"
            gene.is_core = float_frac > core_threshold
            # print("===>", gene)


def extract_genome_id(gene_id):
    """ Extract genome id from gene id. """
    sep = "." if "." in gene_id else ("_" if "_" in gene_id else None)
    if sep is None:
        raise ValueError(f"gene `{gene_id}` does not seem to contain a genome_id.")

    return gene_id[:gene_id.rfind(sep)]


def parse_db_clusters(cluster_data):
    """ Parse gene, cluster, is_core from tsv-stream. """
    return [
        tuple(line.strip().split("\t"))
        for line in get_lines_from_chunks(cluster_data)
    ]


def parse_full_seq_clusters(
    genome_id_prefix,
    genes,
    cluster_data,
    output_dir=None,
):
    """ Parse data from linclust gene clustering. """

    genomes = set()
    cluster_genes = Counter()
    gene_cluster_map = []

    if output_dir is None:
        write_data = True
        cluster_genes_out = nullcontext()
        gene_clusters_out = nullcontext()
    else:
        write_data = False
        cluster_genes_out = open(
            os.path.join(output_dir, f"{genome_id_prefix}.cluster_genes.txt"),
            "wt",
            encoding="UTF-8",
        )
        gene_clusters_out = open(
            os.path.join(output_dir, f"{genome_id_prefix}.gene_clusters.txt"),
            "wt",
            encoding="UTF-8",
        )

    with gene_clusters_out, cluster_genes_out:
        for line in get_lines_from_chunks(cluster_data):
            cluster_id, gene_id = line.split("\t")
            if gene_id.startswith(genome_id_prefix):
                gene_id = gene_id[len(genome_id_prefix) + 1:]

            gene = genes.get(
                gene_id,
                genes.get(gene_id[gene_id.rfind(".") + 1:])
            )
            if gene is not None:
                logger.info("Adding cluster %s to gene %s", cluster_id, gene_id)
                genome_id = gene.genome
                gene_cluster_map.append((gene_id, cluster_id))
            else:
                genome_id = extract_genome_id(gene_id)

            cluster_genes[cluster_id] += 1
            genomes.add(genome_id)
            if write_data:
                print(cluster_id, gene_id, sep="\t", file=cluster_genes_out)
                print(gene_id, cluster_id, sep="\t", file=gene_clusters_out)

        if write_data:
            with open(
                os.path.join(output_dir, f"{genome_id_prefix}.genomes.txt"),
                "wt",
                encoding="UTF-8",
            ) as genomes_out:
                print(*sorted(genomes), sep="\n", file=genomes_out)

    n_genomes = len(genomes)
    return n_genomes, cluster_genes, gene_cluster_map, genome_id_prefix in genomes


@dataclass
class RefGene:
    """ Y-cluster reference gene class. """
    refset: str = None
    speci: int = None
    refset_id: int = None
    is_core: bool = None
    is_singleton: bool = None
    is_unique: bool = None
    prevalence: int = None
    sp100_id: int = None
    genome_id: int = None
    gene_id: str = None
    n_rep_genomes: int = None
    n_rep_genes: int = None

    @classmethod
    def from_string(cls, s):
        """ Construct RefGene from cluster id string. """
        # proMGE095-00037-00002654-AN012-0000748630-000497735_01349-0000000104-0000000608
        fields = s.strip().split("-")
        return cls(
            fields[0],
            int(fields[1]),
            int(fields[2]),
            fields[3][0] == "C",
            fields[3][1] == "S",
            fields[3][1] in "SU",
            int(fields[3][2:]),
            int(fields[4]),
            int(fields[5][:fields[5].find("_")]),
            fields[5],
            int(fields[6]),
            int(fields[7]),
        )


def evaluate_cluster(rep_id, cluster, genes):
    """ Add cluster information from Y-cluster data. """
    ref_genes = {gene for gene in cluster if gene.startswith("proMGE") or gene[0] == "-"}
    query_genes = cluster.difference(ref_genes)

    if query_genes:
        if ref_genes:
            if len(ref_genes) > 1:
                # this is a heuristic --
                # we take the ref_gene with the largest represented genomes ([6])
                # to approximate the cluster's rep_genomes
                ref_genes = sorted(
                    ((int(r.split("-")[6]), r) for r in ref_genes),
                    key=lambda x: x[0], reverse=True
                )
                # print("MULTIREF-CLUSTER", ref_genes[0][1])
                # print(*(r[1] for r in ref_genes[1:]), sep="\n")
                rep = ref_genes[0][1]
            else:
                rep = list(ref_genes)[0]
            ref_gene = RefGene.from_string(rep)
        else:
            # this cluster doesn't contain any reference genes -> genes are accessory
            ref_gene = RefGene(is_core=False, prevalence=0)

        for gene_id in query_genes:
            gene = genes.get(gene_id)
            if gene is not None:
                gene.cluster = rep_id
                gene.is_core = ref_gene.is_core
                gene.prevalence = ref_gene.prevalence


def parse_y_clusters(cluster_data, genes):
    """ Parse data from Y gene clustering approach. """

    cluster, members = None, set()
    for line in get_lines_from_chunks(cluster_data):
        cluster_id, gene_id = line.split("\t")
        if cluster_id != cluster:
            if cluster is not None:
                evaluate_cluster(cluster, members, genes)
            cluster, members = cluster_id, set()

        members.add(gene_id)

    evaluate_cluster(cluster, members, genes)
