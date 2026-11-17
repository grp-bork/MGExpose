# pylint: disable=E0401
""" Module for calculating gene hashes if genes were called externally. """

import pathlib

from ..utils.readers import read_fasta
from ..utils import compute_hash


def hash_genes(args):
	if not args.genes_fasta:
		raise ValueError("Please specify input fasta.")

	compute_gene_hashes(args.genes_fasta, args.genome_id, args.output_dir,)


def compute_gene_hashes(genes_fasta, genome_id, output_dir,):
	outpath = pathlib.Path(output_dir)
	outpath.mkdir(exist_ok=True, parents=True,)
	ffn = outpath / f"{genome_id}.genes.sha256.tsv"
	with open(ffn, "wt", encoding="UTF-8",) as ffn_out:
		for sid, seq in read_fasta(genes_fasta):
			sid = sid[:sid.find(" ")]
			ffn_out.write(f"{sid}\t{compute_hash(seq)}\n")

