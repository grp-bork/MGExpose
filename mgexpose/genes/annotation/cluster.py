import logging

from contextlib import nullcontext

from clustering_parser import parse_db_clusters, parse_full_seq_clusters, parse_y_clusters, evaluate_y_clusters
from ..gene import Gene


logger = logging.getLogger(__name__)


def add_clusters(
	genes,
	fn,
	use_y_clusters=False,
	core_threshold=0.95,
	output_dir=None,
	genome_id=None,
):
	""" Add information from gene clustering to allow for core/accessory gene classification """

	if use_y_clusters:
		if core_threshold == -1:
			parse_y_clusters(fn, genes)
		else:
			evaluate_y_clusters(fn, genes, core_threshold=core_threshold,)
		return None

	write_data = False
	gene_clusters_out = nullcontext()
	n_genomes = 0
	cluster_genes = {}

	with gene_clusters_out:
		if fn is not None:

			if core_threshold != -1:
				n_genomes, cluster_genes, gene_cluster_map, _ = parse_full_seq_clusters(
					genes,
					fn,
					output_dir=output_dir,
				)

				logger.info(
					"Parsed %s genomes with %s gene clusters.",
					n_genomes,
					len(cluster_genes),
				)
			else:
				gene_cluster_map = parse_db_clusters(fn)

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
				if genome_id is None or gene_id.startswith(genome_id):
					gene = genes.get(
						gene_id,
						genes.get(
							gene_id.replace(genome_id + ".", "")
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