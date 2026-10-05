# pylint: disable=E0401
""" Module for genomic island calling. """

import os

from ..genes.annotation.annotate import annotate_genes
from ..islands.island_processing import prepare_islands

from ..utils.writers import dump_islands, write_islands


def island_calling(args):

	annotated_genes, has_clusters = annotate_genes(args)

	if not has_clusters:
		raise ValueError("No clusters found.")

	genomic_islands = prepare_islands(annotated_genes)

	# gff = os.path.join(
	# 	args.output_dir,
	# 	f"{args.genome_id}.genomic_islands.gff3"
	# )

	# dump_islands(genomic_islands, gff)

	gff = os.path.join(args.output_dir, f"{args.genome_id}.genomic_islands.gff3")
	fna = os.path.join(args.output_dir, f"{args.genome_id}.genomic_islands.fna.gz") if args.genome_fasta is not None else None
	write_islands(genomic_islands, gff, genome_fasta=args.genome_fasta, fna=fna,)
