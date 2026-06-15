""" Module docstring """
import os

from ..genes.annotation.annotate import compile_annotations
from ..genes.geneset import GeneSet
from ..islands.annotated_genomic_island import AnnotatedGenomicIsland
from ..islands.genomic_island import GenomicIsland
from ..islands.mge_genomic_island import MgeGenomicIsland
from ..rules.recombinases import get_recombinase_rules
from ..utils.gffio import read_genomic_islands_gff, read_mge_genomic_islands_gff
from ..utils.writers import extract_mge_seqs


def reannotate(args):
    """ docstring """
    if args.annotation_mode == "mges":
        genomic_islands = read_mge_genomic_islands_gff(args.input_genes)
    elif args.annotation_mode == "raw_islands":
        genomic_islands = read_genomic_islands_gff(args.input_genes)
    elif args.annotation_mode == "genes":
        raise NotImplementedError()
    else:
        raise ValueError(f"Unknown annotation_mode: {args.annotation_mode}")

    out_prefix = os.path.join(
        args.output_dir,
        f"{args.genome_id}.mge_islands"
    )

    gene_info_out = open(
        f"{out_prefix}.gene_info.txt", "wt", encoding="UTF-8",
    )

    gene_info_gff = open(
        f"{out_prefix}.gene_info.gff3", "wt", encoding="UTF-8",
    )

    gff_out = open(f"{out_prefix}.gff3", "wt", encoding="UTF-8",)

    mge_islands = {}
    with gene_info_out, gff_out, gene_info_gff:
        print("##gff-version 3", file=gff_out)
        print("##gff-version 3", file=gene_info_gff)

        annotations = compile_annotations(args, multi_run=True,)

        for ct, island in enumerate(genomic_islands):
            # strip previous mge classification from islands (as we want to reannotate)
            island = GenomicIsland.from_island(island, genome_id=args.genome_id,)

            # extract and reannotate the island's genes
            genes = GeneSet.from_island(
                island, composite_gene_ids=args.annotation_mode == "raw_islands",
            )
            _ = genes.annotate(
                annotations, gene_info_out, with_header=ct == 0, gffstream=gene_info_gff,
            )

            # reevaluate annotations and reclassify the island
            island.update_recombinases()
            annotated_island = AnnotatedGenomicIsland.from_island(island)
            mge_island = MgeGenomicIsland.from_island(annotated_island)
            mge_island.evaluate_recombinases(get_recombinase_rules(args.mge_rules))

            mge_island.to_gff(
                gff_out,
                source_db=None,
            )

            mge_islands.setdefault(mge_island.contig, []).append(mge_island)

    if args.genome_fasta:
        extract_mge_seqs(args.genome_fasta, mge_islands, out_prefix)
