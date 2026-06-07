# pylint: disable=E0401
""" Module for gene calling with pyrodigal. """

import pathlib

import pyrodigal

from ..readers import read_fasta


def run_pyrodigal(args):
    """ Call genes with pyrodigal. """
    gf = pyrodigal.GeneFinder(mask=True)

    ids, seqs = zip(*read_fasta(args.genome_fasta))
    _ = gf.train(*seqs)

    outpath = pathlib.Path(args.output_dir)
    outpath.mkdir(exist_ok=True, parents=True,)

    faa_out = open(outpath / f"{args.genome_id}.faa", "wt", encoding="UTF-8",)
    ffn_out = open(outpath / f"{args.genome_id}.ffn", "wt", encoding="UTF-8",)
    gff_out = open(outpath / f"{args.genome_id}.gff", "wt", encoding="UTF-8",)

    with faa_out, ffn_out, gff_out:
        for sid, seq in zip(ids, seqs):
            sid = sid[:sid.find(" ")]
            genes = gf.find_genes(seq)
            genes.write_translations(faa_out, sid)
            genes.write_genes(ffn_out, sid)
            genes.write_gff(gff_out, sid, full_id=False,)
