import pathlib

import pyrodigal

from .readers import read_fasta


def run_pyrodigal(args):
	gf = pyrodigal.GeneFinder(mask=True)

	ids, seqs = zip(*read_fasta(args.genome_fasta))
	_ = gf.train(*seqs)

	outpath = pathlib.Path(args.output_dir)
	outpath.mkdir(exist_ok=True, parents=True,)

	faa_out = open(outpath / f"{args.genome_id}.faa", "wt")
	ffn_out = open(outpath / f"{args.genome_id}.ffn", "wt")
	gff_out = open(outpath / f"{args.genome_id}.gff", "wt")

	with faa_out, ffn_out, gff_out:
		for i, (sid, seq) in enumerate(zip(ids, seqs), start=1):
			sid = sid[:sid.find(" ")]
			genes = gf.find_genes(seq)
			genes.write_translations(faa_out, sid)
			genes.write_genes(ffn_out, sid)
			genes.write_gff(gff_out, sid, full_id=False,)

	# #with pyhmmer.easel.SequenceFile("/g/bork7/fullam/progenomes/pg4/download/bacteria/ncbi_dataset/data/GCA_000005825.2/GCA_000005825.2_ASM582v2_genomic.fna.gz") as seq_file:
	# with pyhmmer.easel.SequenceFile("/g/bork7/fullam/progenomes/pg4/download/bacteria/ncbi_dataset/data/GCA_942641705.2/GCA_942641705.2_Pse-ATUE_S32H133_draftAssembly_v2_genomic.fna.gz") as seq_file:
	# 	seqs = list(seq_file)


	# training_info = gf.train(*seqs)
	# print(training_info)
	# print(seqs)
	# print(seqs[0])


	# #genes = gf.find_genes(*seqs)
	# #print(genes)

	# with open("test.ffn", "wt") as ffn_out, open("test.faa", "wt") as faa_out, open("test.gff", "wt") as gff_out:
	# 	for i,seq in enumerate(seqs,start=1):
	# 		genes = gf.find_genes(seq)
	# 		genes.write_translations(faa_out, f"XYZ.{i}")#, "XYZ")
	# 		genes.write_genes(ffn_out, f"XYZ.{i}")#, "XYZ")
	# 		genes.write_gff(gff_out, f"XYZ.{i}", full_id=False,)#, "XYZ")

	# #for gene in genes:
	# #    print(gene)
	# #    break
