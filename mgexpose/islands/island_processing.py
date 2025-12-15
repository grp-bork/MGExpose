import logging

from .annotated_genomic_island import AnnotatedGenomicIsland
from .genomic_island import GenomicIsland
from .mge_genomic_island import MgeGenomicIsland


logger = logging.getLogger(__name__)


def compute_genomic_islands(genes):
    contigs = {}
    for gene in genes:
        if gene.has_basic_annotation():
            contigs.setdefault((gene.speci, gene.contig), []).append(gene)

    for _, genes in sorted(contigs.items()):
        current_island = None

        for gene in sorted(genes, key=lambda g: (g.start, g.end, g.strand)):
            boundary_found = (
                current_island is None or 
                not current_island.has_compatible_genome_type(gene)
            )
            if boundary_found:
                if current_island is not None:
                    if current_island.recombinases:
                        yield current_island

                current_island = GenomicIsland.from_gene(gene)

            current_island.add_gene(gene)

        if current_island is not None:
            if current_island.recombinases:
                logger.info("GenomicIsland %s created.", str(current_island))
                yield current_island


def add_genes_to_precomputed_islands(genes, islands):
    logger.info("Precomputed islands: %s", len(islands))
    for gene in genes:
        is_annotated = gene.has_basic_annotation(skip_core_gene_computation=True)
        add_gene = False
        if is_annotated:
            gene.contig = gene.contig.split(".")[-1]
            for island in islands.get(gene.contig, []):
                log_str = (
                    f"Checking gene={gene.contig}:"
                    f"{gene.start}-{gene.end} against "
                    f"{island.contig}:{island.start}-{island.end}: "
                )

                if gene.is_in_interval(island.start, island.end):
                    add_gene = True
                    if gene.speci is None or gene.speci == "no_speci":
                        gene.speci = {island.name}
                    else:
                        gene.speci.add(island.name)
                    gene.parent = island.get_id()
                    island.add_gene(gene)

                logger.info("%s %s", log_str, str(add_gene))
    
    for _, islands in islands.items():
        for island in islands:
            if island.recombinases:
                logger.info("GenomicIsland %s created.", str(island))
                yield island
        

def prepare_islands(genes, precomputed_islands=None,):

    if precomputed_islands is None:
        islands = compute_genomic_islands(genes)
    else:
        islands = add_genes_to_precomputed_islands(genes, precomputed_islands)

    yield from islands


def annotate_islands(islands,):
    """ Adds annotation to previously computed islands. """
    for island in sorted(islands, key=lambda x: (x.contig, x.start, x.end)):
        annotated_island = AnnotatedGenomicIsland.from_genomic_island(island)
        yield annotated_island


def evaluate_islands(islands, rules,):
    """ Classify/annotate mge islands according to present signals. """
    for island in islands:
        mge_island = MgeGenomicIsland.from_annotated_genomic_island(island)
        mge_island.evaluate_recombinases(rules)
        yield mge_island