# pylint: disable=E0401
""" Module for gene calling with pyrodigal. """

import pathlib

import pyrodigal

from ..utils.readers import read_fasta


def gene_calling(args):

    if not args.genome_fasta:
        raise ValueError("Please specify genome input.")

    run_pyrodigal(args.genome_fasta, args.genome_id, args.output_dir, pr_meta=args.meta)

def run_pyrodigal(genome_fasta, genome_id, output_dir, pr_meta=False,):
    """ Call genes with pyrodigal. """
    gf = pyrodigal.GeneFinder(mask=True, meta=pr_meta,)

    ids, seqs = zip(*read_fasta(genome_fasta))
    _ = gf.train(*seqs)

    outpath = pathlib.Path(output_dir)
    outpath.mkdir(exist_ok=True, parents=True,)

    faa = outpath / f"{genome_id}.faa"
    ffn = outpath / f"{genome_id}.ffn"
    gff = outpath / f"{genome_id}.gff"

    faa_out = open(faa, "wt", encoding="UTF-8",)
    ffn_out = open(ffn, "wt", encoding="UTF-8",)
    gff_out = open(gff, "wt", encoding="UTF-8",)

    with faa_out, ffn_out, gff_out:
        for sid, seq in zip(ids, seqs):
            sid = sid[:sid.find(" ")]
            genes = gf.find_genes(seq)
            genes.write_translations(faa_out, sid)
            genes.write_genes(ffn_out, sid)
            genes.write_gff(gff_out, sid, full_id=False,)

    return faa, ffn, gff
