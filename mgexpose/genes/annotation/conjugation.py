""" Module to add conjugation system annotations. """


def add_conjugation_systems(conjugation_systems, genes,):
    """ Add information from txsscan """
    for gene_id, conjugation_data in conjugation_systems:
        gene = genes.get(gene_id)
        if gene is not None:
            for sgene, system, rule, *_ in conjugation_data:
                gene.conjugation_systems.append(f"{sgene}:{system}")
                gene.conjugation_rules.append(rule)
