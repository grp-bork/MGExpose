import os

from ..genes.annotation.annotate import compile_annotations
from ..genes.geneset import GeneSet
from ..islands.island_processing import annotate_islands, evaluate_islands, prepare_islands
from ..rules import get_recombinase_rules
from ..utils.readers import read_precomputed_islands
from ..utils.writers import extract_mge_seqs, write_final_results


def denovo(args):
	genomic_islands = None

	genes = GeneSet.from_file(
		args.input_genes,
		genome_id=args.genome_id,
		species=args.species,
		gene_type=args.input_gene_type,
		composite_gene_ids=args.dbformat != "PG3",
	)

	annotations, has_clusters = compile_annotations(args, multi_run=False,)

	gene_info_out = open(
		os.path.join(
			args.output_dir,
			f"{args.genome_id}.gene_info.txt",
		),
		"wt",
		encoding="UTF-8",
	)

	gene_info_gff = open(
		os.path.join(
			args.output_dir,
			f"{args.genome_id}.gene_info.gff3",
		),
		"wt",
		encoding="UTF-8",
	)


	with gene_info_out, gene_info_gff:
		annotated_genes = genes.annotate(annotations, stream=gene_info_out, gffstream=gene_info_gff,)

		precomputed_islands = None
		if args.single_island or args.precomputed_islands:
			precomputed_islands = read_precomputed_islands(
				genome_id=args.genome_id,
				single_island=args.single_island,
				island_file=args.precomputed_islands,
			)

		if has_clusters or precomputed_islands or args.contigs_are_islands:
			genomic_islands = prepare_islands(
				annotated_genes, precomputed_islands=precomputed_islands, contigs_are_islands=args.contigs_are_islands,
			)

			annotated_islands = annotate_islands(genomic_islands,)

			mge_islands = list(
				evaluate_islands(annotated_islands, get_recombinase_rules(args.mge_rules),)
			)

			if mge_islands:
				write_final_results(mge_islands, args)