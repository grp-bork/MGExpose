
""" Module docstring """
import os
import pathlib

from ..genes.geneset import GeneSet
from ..islands.annotated_genomic_island import AnnotatedGenomicIsland
from ..islands.genomic_island import GenomicIsland
from ..islands.mge_genomic_island import MgeGenomicIsland
from ..rules.recombinases import get_recombinase_rules
from ..utils.gffio import read_mge_genomic_islands_gff
from ..utils.writers import extract_mge_seqs


def liftover(args):
	""" docstring """
	mge_islands = {}
	mge_rules = get_recombinase_rules(args.mge_rules)
	# with open(args.island_mapping, "rt") as _in:
	#     # island mapping format:
	#     # source -> dest, hence we need to reverse
	#     # as we're checking the dest islands only
	#     island_mapping = dict(
	#         line.strip().split()[::-1]
	#         for line in _in
	#     )
	if args.island_mapping is None:
		raise ValueError("No islands to map specified.")
	island_mapping = dict()
	if pathlib.Path(args.island_mapping).is_file():
		with open(args.island_mapping, "rt") as _in:
			island_mapping = dict(
				l.strip().split(",")[::-1]
				for l in _in
			)
	elif args.island_mapping != "all":
		island_mapping = dict([args.island_mapping.strip().split(",")[::-1]])
	print("ISLAND_MAPPING", island_mapping)
	source_islands = {
		island.get_id(): island
		for island in read_mge_genomic_islands_gff(args.source_islands)
	}
	dest_islands = {
		island.get_id(): island
		for island in read_mge_genomic_islands_gff(args.dest_islands)
	}
	out_prefix = os.path.join(
		args.output_dir,
		f"{args.genome_id}.mge_islands.liftover"
	)

	out_gff3 = f"{out_prefix}.gff3"
	i = 1
	while os.path.isfile(out_gff3):
		out_gff3 = f"{out_gff3}.{i}"
		i += 1

	gff_out = open(out_gff3, "wt", encoding="UTF-8",)

	print("source_islands", *source_islands, sep="\n")
	print("dest_islands", *dest_islands, sep="\n")

	with gff_out:
		print("##gff-version 3", file=gff_out)
		for id1, dst in dest_islands.items():
			id2 = island_mapping.get(id1) if island_mapping else id1
			src = source_islands.get(id2)
			print(f"{id1=}, {id2=}, {src=}")
			new_island = GenomicIsland.from_island(dst, dst.genome,)
			if src is None:
				new_island.update_recombinases()
				annotated_island = AnnotatedGenomicIsland.from_island(new_island)
				mge_island = MgeGenomicIsland.from_island(annotated_island)
				mge_island.evaluate_recombinases(mge_rules)
			else:
				# new_island = GenomicIsland.from_island(dst, dst.genome,)
				# for src_gene, dst_gene in zip(src.genes, new_island.genes):
				#     dst_gene.liftover(src_gene)
				GeneSet.liftover(tuple(src.get_genes()), tuple(new_island.get_genes()))
				new_island.update_recombinases()
				annotated_island = AnnotatedGenomicIsland.from_island(new_island)
				mge_island = MgeGenomicIsland.from_island(annotated_island)
				mge_island.evaluate_recombinases(mge_rules)
			mge_island.to_gff(
				gff_out,
				source_db=None,
			)
			mge_islands.setdefault(mge_island.contig, []).append(mge_island)

	if args.extract_islands:
		extract_mge_seqs(args.extract_islands, mge_islands, out_prefix)