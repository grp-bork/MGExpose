#!/usr/bin/env python
# pylint: disable=R0912,R0914,R0915
import os

from collections import Counter, defaultdict

from .gffio import read_mge_genomic_islands_gff #  FIXME

from .base_logger import logger


''' Quick and dirty implementation - actually does not belong here ''' # FIXME
def stat_core(gff_files):
    """Calculate statistics on MGEs in the core genome and return as a dictionary."""
    total = 0
    core = 0
    for f in gff_files:
        island = read_mge_genomic_islands_gff(f)
        if island.is_core:
            core += 1
        total += 1
    return {
        "total": total,
        "core_count": core,
        "accessory_count": total - core,
    }

def stat_nested(islands):
    """Calculate statistics on nested MGEs and return as a dictionary."""
    total = 0
    nested = 0
    for island in islands:
        if island.mge_type == "nested":
            nested += 1
        total += 1
    nested_percentage = (nested / total * 100) if total > 0 else 0
    non_nested_percentage = ((total - nested) / total * 100) if total > 0 else 0

    return {
        "total": total,
        "nested_count": nested,
        "non_nested_count": total - nested,
        "nested_percentage": nested_percentage,
        "non_nested_percentage": non_nested_percentage,
    }


def count_nested(islands):
    """
    Calculate the count of nested MGEs.
    
    Parameters:
    - islands: List of MGE objects.
    
    Returns:
    - Integer count of nested MGEs.
    """
    nested = sum(1 for island in islands if island.mge_type == "nested")
    return nested


def count_core(islands):
    """
    Calculate the count of MGEs in the core genome.
    
    Parameters:
    - islands: List of MGE objects.
    
    Returns:
    - Integer count of MGEs in the core genome.
    """
    core = sum(1 for island in islands if island.is_core)
    return core


def count_total_islands(islands):
    """
    Count the total number of MGE islands.
    
    Parameters:
    - islands: List of MGE objects.
    
    Returns:
    - Integer count of total MGE islands.
    """
    return len(list(islands))


def stat_mge_type(islands):
    """
    Calculate counts of each MGE type and return as a Counter object.
    
    Parameters:
    - islands: List of MGE objects.
    
    Returns:
    - Counter object with counts of each MGE type.
    """
    mge_counts = Counter()
    for island in islands:
        try:
            if island.mge_type == "nested":
                mge_counts["nested"] += 1
            else:
                # Get the first key (assuming there's only one key)
                mge = next(iter(island.mge.keys()))
                mge_counts[mge] += 1
        except Exception as e:
            raise ValueError(f"Unknown or absent MGE type: {e}")
    
    return mge_counts


def stat_mean_genes(islands):
    """Calculate the mean number of genes per MGE and return as a dictionary."""
    genes_lst = [island.n_genes for island in islands]
    mean_genes = (sum(genes_lst) / len(genes_lst)) if genes_lst else 0

    # Return the result as a dictionary
    return {"mean_genes_per_mge": mean_genes}


def extract_cargo(island):
    cargo_genes = []
    for gene in island.genes:
        if (gene.phage is None) and (gene.recombinase is None) and (gene.secretion_system is None):
            cargo_genes.append(gene)
    return cargo_genes


def get_kegg_ko(gene):
    for key, value in gene.eggnog:
        if key == "kegg_ko":
            return value


def get_cazy(gene):
    for key, value in gene.eggnog:
        if key == "cazy":
            return value


def get_cog_category(gene):
    if gene.eggnog:
        for key, value in gene.eggnog:
            if key == "cog_fcat":
                if value:
                    return value
                else:
                    return '-'
            else:
                return 'no_cog_fcat'
    else:
        return 'no_eggnog'


# SP95 for SPIRE 
def get_gene_cluster(gene):
    return gene.cluster


def get_gene_id(gene):
    id_lst = gene.id.split('.') # orginal ID e.g. 28901.SAMN15849311.GCA_014242155_04079
    return id_lst[2] # extract to match CY clustering IDs e.g. GCA_xxx_xxx     

# Extract MGE recombinases
def extract_mger(island):
    mgeR_genes = []
    for gene in island.genes:
        if gene.recombinase:
            try:
                gene_id = gene.id
                gene_cluster = get_gene_cluster(gene)
                mgeR = gene.recombinase 
                annot = [gene_id, gene_cluster, mgeR]
                mgeR_genes.append(annot)
            except Exception as e:
                logger.error(f"Error processing recombinase gene {gene}: {e}")
                logger.error(traceback.format_exc())
    return mgeR_genes


# Extract secretion system genes
def extract_secretion_system(island):
    secretion_genes = []
    for gene in island.genes:
        if gene.secretion_system:
            try:
                gene_id = gene.id
                gene_cluster = get_gene_cluster(gene)
                info = gene.secretion_system
                annot = [gene_id, gene_cluster, info]
                secretion_genes.append(annot)
            except Exception as e:
                logger.error(f"Error processing secretion system gene {gene}: {e}")
                logger.error(traceback.format_exc())
    return secretion_genes


# Extract phage genes
def extract_phage(island):
    phage_genes = []
    for gene in island.genes:
        if gene.phage:
            try:
                gene_id = gene.id
                gene_cluster = get_gene_cluster(gene)
                info = gene.phage
                annot = [gene_id, gene_cluster, info]
                phage_genes.append(annot)
            except Exception as e:
                logger.error(f"Error processing phage gene {gene}: {e}")
                logger.error(traceback.format_exc())
    return phage_genes


def get_most_common_kegg_ko(genes):
    """
    Calculate the most common KEGG KO annotations for a list of genes.
    
    Parameters:
    - genes: List of gene objects.
    
    Returns:
    - Counter object with counts of KEGG KO annotations.
    """
    kos = [get_kegg_ko(gene) for gene in genes]
    return Counter(kos)


def count_gene_clusters(genes, **kwargs):
    """
    Count genes in gene clusters.
    
    Parameters:
    - genes: List of gene objects.
    
    Returns:
    - Counter object with counts of cluster genes.
    """
    gene_clusters = [get_gene_cluster(gene) for gene in genes]
    return Counter(gene_clusters)


def get_majority_cog_category(genes, **kwargs):
    cluster_to_categories = defaultdict(list)
    
    for gene in genes:
        gene_cluster = get_gene_cluster(gene)
        cog_category = get_cog_category(gene)
        cluster_to_categories[gene_cluster].append(cog_category)
    
    majority_cog_category = {}
    
    # Determine the majority cog category for each cluster
    for cluster, categories in cluster_to_categories.items():
        category_count = Counter(categories)
        majority_cog_category[cluster] = category_count.most_common(1)[0][0] # a list of tuples where each tuple is a category and its count; [0][0] extracts the category with the highest count.
    #logger.info(majority_cog_category)
    return majority_cog_category


def get_gene_annotation(genes, func, **kwargs):
    gene_annot_dict = {}
    for gene in genes:
        gene_id = get_gene_id(gene)
        func_annot = func(gene)
        gene_annot_dict[gene_id] = func_annot
        
    return gene_annot_dict


def get_genes_cog_categories(genes, **kwargs):
    return get_gene_annotation(genes, get_cog_category)


def get_genes_ko(genes, *_args):
    return get_gene_annotation(genes, get_kegg_ko)


def get_genes_cazy(genes, *_args):
    return get_gene_annotation(genes, get_cazy)


def get_genes_clusters(genes, **kwargs):
    return get_gene_annotation(genes, get_gene_cluster)


def get_genes(genes, mge_id):
    gene_ids = [get_gene_id(gene) for gene in genes]
    return {mge_id: gene_ids}
    

def count_per_mge_cargo(islands, func):
    """
    Extract and count cargo genes associated with each MGE type and return as a dictionary of Counter objects.
    
    Parameters:
    - islands: List of MGE objects.
    - func: function to extract cargo e.g. count KO terms or gene clusters 
    
    Returns:
    - Dictionary where keys are MGE types, and values are Counter objects with KEGG KO counts.
    """
    mge_cargo_counts = {
        "nested": Counter(),
        "phage": Counter(),
        "phage_like": Counter(),
        "is_tn": Counter(),
        "ce": Counter(),
        "mi": Counter(),
        "integron": Counter(),
        "cellular": Counter(),
    }

    for island in islands:
        try:
            cargo = extract_cargo(island)
            if island.mge_type == "nested":
                mge_cargo_counts["nested"].update(func(cargo))
            else:
                # Get the first key (assuming there's only one key)
                mge = next(iter(island.mge.keys()))
                mge_cargo_counts[mge].update(func(cargo))
        except Exception as e:
            raise ValueError(f"Error processing cargo for island: {e}")
    
    return mge_cargo_counts

    
def get_per_mge_cargo(islands, func):
    mge_cargo_annot = {
        "nested": {},
        "phage": {},
        "phage_like": {},
        "is_tn": {},
        "ce": {},
        "mi": {},
        "integron": {},
        "cellular": {},
    }

    for island in islands:
        try:
            cargo = extract_cargo(island)
            mge_id = island.get_id()
            if island.mge_type == "nested":
                mge_cargo_annot["nested"].update(func(cargo, mge_id)) # Merges another dictionary into existing one
            else:
                # Get the first key (assuming there's only one key)
                mge = next(iter(island.mge.keys()))
                mge_cargo_annot[mge].update(func(cargo, mge_id))
        except Exception as e:
            raise ValueError(f"Error processing cargo for island: {e}")
    
    return mge_cargo_annot


def get_machinery_genes_tsv(islands):
    tsv_rows = []
    header = ['mge_id', 'mge', 'n_genes', 'gene_id', 'gene_cluster', 'feature_type', 'feature_info']
    tsv_rows.append('\t'.join(header))

    for island in islands:
        try:
            mge_id = island.get_id()
            n_genes = len(island.genes)
            mge = ",".join(f"{k}:{v}" for k, v in island.mge.items())

            recombinases = extract_mger(island)  # list of [gene_id, gene_cluster, info]
            conj_machinery = extract_secretion_system(island)  
            phage_machinery = extract_phage(island)

            for gene in recombinases:
                gene_id, gene_cluster, info = gene
                tsv_rows.append(f"{mge_id}\t{mge}\t{n_genes}\t{gene_id}\t{gene_cluster}\tmgeR\t{info}")

            for gene in conj_machinery:
                gene_id, gene_cluster, info = gene
                tsv_rows.append(f"{mge_id}\t{mge}\t{n_genes}\t{gene_id}\t{gene_cluster}\tsecretion_system\t{info}")

            for gene in phage_machinery:
                gene_id, gene_cluster, info = gene
                tsv_rows.append(f"{mge_id}\t{mge}\t{n_genes}\t{gene_id}\t{gene_cluster}\tphage\t{info}")

        except Exception as e:
            raise ValueError(f"Error processing machinery for island {island}: {e}")

    tsv_output = '\n'.join(tsv_rows)
    return tsv_output

def get_cargo_genes_tsv(islands):
    tsv_rows = []
    header = ['mge_id', 'mge', 'gene_ids', 'gene_clusters']
    tsv_rows.append('\t'.join(header))

    for island in islands:
        gene_ids = []
        gene_clusters = []
        try:
            mge_id = island.get_id()
            mge = ",".join(f"{k}:{v}" for k, v in island.mge.items())

            cargo_genes = extract_cargo(island)

            for gene in cargo_genes:
                gene_ids.append(gene.id)
                gene_clusters.append(get_gene_cluster(gene))
            gene_ids = ';'.join(gene_ids)
            gene_clusters = ';'.join(gene_clusters)
            tsv_rows.append(f"{mge_id}\t{mge}\t{gene_ids}\t{gene_clusters}")

        except Exception as e:
            raise ValueError(f"Error processing cargo for island {island}: {e}")

    tsv_output = '\n'.join(tsv_rows)
    return tsv_output

# Counting works with aggregation since the objects are small. Getting only works with batch saving, since the output is huge.
def count_cargo_gene_clusters(islands):
    return count_per_mge_cargo(islands, count_gene_clusters)

    
def get_cargo_species_gene_clusters(islands):
    return get_per_mge_cargo(islands, func=get_genes_clusters)


# Output: for each mge_type output a dictionary with geneID: COG. geneID supposed to be unique -> overwriting is okay.
def get_cargo_genes_cog(islands):
    return get_per_mge_cargo(islands, func=get_genes_cog_categories)


def get_cargo_genes_ko(islands):
    return get_per_mge_cargo(islands, func=get_genes_ko)


def get_cargo_genes_cazy(islands):
    return get_per_mge_cargo(islands, func=get_genes_cazy)


# Output: for each mge_type output a dictionary mge_id: list(cargo_ids)
def get_cargo_genes(islands):
    return get_per_mge_cargo(islands, func=get_genes)
    
    
