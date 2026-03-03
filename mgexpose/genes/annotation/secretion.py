""" Module to add secretion system annotations. """


def add_secretion_systems(secretion_systems, genes,):
    """ Add information from txsscan """
    for gene_id, secretion_data in secretion_systems:
        gene = genes.get(gene_id)
        if gene is not None:
            for sgene, system, rule, *_ in secretion_data:
                gene.secretion_systems.append(f"{sgene}:{system}")
                gene.secretion_rules.append(rule)
