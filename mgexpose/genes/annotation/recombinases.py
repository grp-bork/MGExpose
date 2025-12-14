from ..readers import read_recombinase_hits


def add_recombinases(genes, fn, pyhmmer=False,):
	""" Add information from recombinase scan """
	for gene_id, recombinase in read_recombinase_hits(fn, pyhmmer=pyhmmer,):
		gene = genes.get(gene_id)
		if gene is not None:
			gene.recombinase = recombinase