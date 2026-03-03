""" Module to add recombinase annotations. """


def add_recombinases(recombinases, genes,):
    """ Add information from recombinase scan """
    for gene_id, recombinase in recombinases:
        print("RECOMBINASE", gene_id, recombinase)
        gene = genes.get(gene_id)
        if gene is not None:
            gene.recombinase = recombinase
