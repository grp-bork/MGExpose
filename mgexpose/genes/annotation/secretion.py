""" Module to add secretion system annotations. """
from ...utils.readers import parse_macsyfinder_report


def add_secretion_systems(fn, genes, rules_fn):
    """ Add information from txsscan """
    for gene_id, secretion_data in parse_macsyfinder_report(fn, rules_fn):
        gene = genes.get(gene_id)
        if gene is not None:
            for sgene, system, rule, *_ in secretion_data:
                gene.secretion_systems.append(f"{sgene}:{system}")
                gene.secretion_rules.append(rule)
