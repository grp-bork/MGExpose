""" Module for processing mobile genetic islands """

import contextlib
import logging


from .islands import GenomicIsland, AnnotatedGenomicIsland, MgeGenomicIsland


logger = logging.getLogger(__name__)


def is_valid_stream(stream):
    """ Checks if a stream-variable represents a valid stream. """
    return stream is not None and not isinstance(
        stream, contextlib.nullcontext
    )


def check_island_genes(genes, precomputed_islands=None):
    """ Check if genes have valid annotation and belong to a precomputed island. """
    has_precomputed_islands = False
    if precomputed_islands is not None:
        logger.info("Precomputed islands: %s", len(precomputed_islands))
        has_precomputed_islands = True

    for gene in genes:
        is_annotated = gene.has_basic_annotation(skip_core_gene_computation=has_precomputed_islands)
        add_gene = False
        if is_annotated and has_precomputed_islands:
            gene.contig = gene.contig.split(".")[-1]
            for island in precomputed_islands.get(gene.contig, []):
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
        else:
            add_gene = is_annotated

        log_str = (
            f"Gene {gene}: {is_annotated=} {add_gene=} ->"
            f" contig set = {is_annotated and add_gene}"
        )
        logger.info(log_str)
        if add_gene:
            yield gene


def filter_precomputed_islands(precomputed_islands, raw_island_stream=None, island_stream=None):
    """ Return precomputed islands with mge signals. """
    for _, islands in precomputed_islands.items():
        seen_islands = set()
        for island in islands:
            island.dump(seen_islands, raw_outstream=raw_island_stream, outstream=island_stream)
            if island.recombinases:
                logger.info("GenomicIsland %s created.", str(island))
                yield island


def compute_islands(contigs, raw_island_stream=None, island_stream=None):
    """ Form genomic islands from stretches of core or accessory genes,
    the return those with mge signals."""
    for _, genes in sorted(contigs.items()):
        seen_islands = set()
        current_island = None
        for gene in sorted(genes, key=lambda g: (g.start, g.end, g.strand)):
            if current_island is None or current_island.is_core != gene.is_core:
                if current_island is not None:
                    current_island.dump(
                        seen_islands,
                        raw_outstream=raw_island_stream,
                        outstream=island_stream
                    )
                    if current_island.recombinases:
                        yield current_island

                current_island = GenomicIsland.from_gene(gene)

            current_island.add_gene(gene)

        if current_island is not None:
            current_island.dump(
                seen_islands,
                raw_outstream=raw_island_stream,
                outstream=island_stream
            )
            if current_island.recombinases:
                yield current_island


def generate_island_set(
    genes,
    pang_calls_out=None,
    raw_islands_out=None,
    islands_out=None,
    precomputed_islands=None,
):
    """ Compute mge islands """
    contigs = {}
    logger.info("generate_island_set: collecting genes")

    for gene in check_island_genes(genes, precomputed_islands=precomputed_islands):
        logger.info("Adding gene %s to contig set.", str(gene))
        if is_valid_stream(pang_calls_out):
            print(gene, file=pang_calls_out)
        if precomputed_islands is None:
            contigs.setdefault(
                (gene.speci, gene.contig), []
            ).append(gene)

    if is_valid_stream(islands_out):
        print(*GenomicIsland.get_fieldnames(), sep="\t", file=islands_out)
        island_stream = islands_out
    else:
        island_stream = islands_out
    if is_valid_stream(raw_islands_out):
        print(*GenomicIsland.RAW_TABLE_HEADER, sep="\t", file=raw_islands_out)
        raw_island_stream = raw_islands_out
    else:
        raw_island_stream = None

    if precomputed_islands is not None:
        yield from filter_precomputed_islands(
            precomputed_islands,
            raw_island_stream=raw_island_stream,
            island_stream=island_stream,
        )
    else:
        yield from compute_islands(
            contigs,
            raw_island_stream=raw_island_stream,
            island_stream=island_stream,
        )


def annotate_islands(islands, outstream=None):
    """ Adds annotation to previously computed islands. """
    do_print = is_valid_stream(outstream)
    if do_print:
        print(*AnnotatedGenomicIsland.TABLE_HEADERS, sep="\t", file=outstream)
    for island in sorted(islands, key=lambda x: (x.contig, x.start, x.end)):
        annotated_island = AnnotatedGenomicIsland.from_genomic_island(island)
        if do_print:
            print(annotated_island, file=outstream)
        yield annotated_island


def evaluate_islands(islands, rules, outstream=None, outstream2=None):
    """ Classify/annotate mge islands according to present signals. """
    if is_valid_stream(outstream):
        print(*MgeGenomicIsland.TABLE_HEADERS, sep="\t", file=outstream)
    if is_valid_stream(outstream2):
        print(*MgeGenomicIsland.TABLE_HEADERS, sep="\t", file=outstream2)
    for island in islands:
        mge_island = MgeGenomicIsland.from_annotated_genomic_island(island)
        mge_island.evaluate_recombinases(
            rules,
            outstream=outstream if not isinstance(outstream, contextlib.nullcontext) else None,
            outstream2=outstream2 if not isinstance(outstream2, contextlib.nullcontext) else None,
        )

        yield mge_island


def prepare_precomputed_islands(single_island=None, island_file=None, genome_id=None,):
    """ Helper function to deal with precomputed regions/islands. """
    precomputed_islands = None
    if single_island and island_file:
        raise ValueError("Both --single_island and --precomputed_islands set.")
    if single_island and not island_file:
        precomputed_islands = [GenomicIsland.from_region_string(single_island, genome_id=genome_id,)]
    elif not single_island and island_file:
        with open(island_file, "rt", encoding="UTF-8",) as _in:
            precomputed_islands = [GenomicIsland.from_region_string(line, genome_id=genome_id,) for line in _in]

    if precomputed_islands is not None:
        precomputed_islands_by_contig = {}
        for island in precomputed_islands:
            precomputed_islands_by_contig.setdefault(island.contig, []).append(island)
        precomputed_islands = precomputed_islands_by_contig

    return precomputed_islands
