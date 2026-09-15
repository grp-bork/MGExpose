# pylint: disable=E0401
""" Module for gene calling with pyrodigal. """

import hashlib
import pathlib

from io import StringIO

import pyrodigal

from ..utils.readers import read_fasta

def reverse_complement(sequence):
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
    }
    return "".join(complement[base] for base in sequence.upper()[::-1])


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
        has_header = False
        for sid, seq in zip(ids, seqs):
            sid = sid[:sid.find(" ")]
            genes = gf.find_genes(seq)
            genes.write_translations(faa_out, sid)
            genes.write_genes(ffn_out, sid)

            buf = StringIO()
            
            genes.write_gff(buf, sid, full_id=False,)

            # genes.write_gff(gff_out, sid, full_id=False,)
            # buf.seek(len(buf.getvalue()) - 1)
            gfflines = buf.getvalue().split("\n")
            if not has_header:
                gff_out.write(f"{gfflines[0]}\n")
                has_header = True
            gff_out.write(f"{gfflines[1]}\n")
            gff_out.write(f"{gfflines[2]}\n")
            for gene, line in zip(genes, gfflines[3:]):
                fwd = hashlib.sha256(gene.sequence().encode()).hexdigest()
                rev = hashlib.sha256(reverse_complement(gene.sequence()).encode()).hexdigest()
                gff_out.write(f"{line[:-1]};fwd={fwd};rev={rev}\n")
                
                



    return faa, ffn, gff
