""" Module to add recombinase annotations. """
from ...utils.readers import read_recombinase_hits


def add_recombinases(fn, genes, pyhmmer=False,):
    """ Add information from recombinase scan """
    for gene_id, recombinase in read_recombinase_hits(fn, pyhmmer=pyhmmer,):
        print("RECOMBINASE", gene_id, recombinase)
        gene = genes.get(gene_id)
        if gene is not None:
            gene.recombinase = recombinase
