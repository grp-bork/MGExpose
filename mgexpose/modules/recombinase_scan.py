# pylint: disable=C0301,E0401,R0914
""" Module to detect and annotate recombinases via pyhmmer. """

import pathlib

import pyhmmer

from .gene_calling import run_pyrodigal
from ..recombinases import MGE_ALIASES
from ..rules.recombinases import get_recombinase_rules
from ..utils.gffio import read_prodigal_gff


RECOMBINASE_SCAN_HEADER = (
    "#unigene",
    "recombinase_SMART_hmm_name",
    "PFAM_accession",
    "MGE_prediction",
    "hmmsearch_fullsequence_evalue",
    "hmmsearch_fullsequence_score",
    "MGE_prediction_confidence",
)


def get_protein_coords(gff):
    """ Read protein coordinates from gff. """

    proteins = {}
    for gene in read_prodigal_gff(gff):
        gene.id = f'{gene.contig}_{gene.id.split("_")[-1]}'
        proteins[gene.id] = gene
    return proteins



def recombinase_scan(args):
    faa, gff = args.protein_fasta, args.protein_coords
    has_genome_input = bool(args.genome_fasta)
    has_proteome_input = faa and gff

    if has_genome_input and has_proteome_input:
        raise ValueError("Please specify either proteome input or genome input. Not both.")
    if not has_genome_input and not has_proteome_input:
        raise ValueError("Please specify either genome input (--genome_fasta) or proteome input (--protein_fasta & --protein_coords).")

    if has_genome_input:
        faa, _, gff = run_pyrodigal(args.genome_fasta, args.genome_id, args.output_dir, pr_meta=args.meta,)
    
    run_pyhmmer(faa, gff, args.genome_id, args.recombinase_hmms, args.output_dir,
                mge_rules=args.mge_rules, threads=args.threads,)


def run_pyhmmer(faa, gff, genome_id, recombinase_hmms, output_dir, mge_rules=None, threads=1,):
    """ Detect and annotate recombinases via pyhmmer. """

    proteins = get_protein_coords(gff)

    mge_rules = get_recombinase_rules(mge_rules, for_recombinase_scan=True,)

    with pyhmmer.easel.SequenceFile(
        faa,
        digital=True,
        alphabet=pyhmmer.easel.Alphabet.amino(),
    ) as seq_file:
        protein_seqs = list(seq_file)
    with pyhmmer.plan7.HMMFile(recombinase_hmms) as hmm_file:
        hmm_hits = list(
            pyhmmer.hmmsearch(
                hmm_file,
                protein_seqs,
                cpus=threads,
                backend="multiprocessing",
                bit_cutoffs="gathering",
            )
        )

    outpath = pathlib.Path(output_dir)
    outpath.mkdir(exist_ok=True, parents=True,)

    raw_table_out = open(
        outpath / f"{genome_id}.recombinase_hmmsearch.out",
        "wb"
    )

    with raw_table_out:
        seen = {}
        for i, hits in enumerate(hmm_hits):
            write_header = i == 0
            hits.write(raw_table_out, header=write_header)
            for hit in hits:
                hit_name = hit.name
                for domain in hit.domains:
                    best_score = seen.setdefault(hit_name, (0.0, None, None))[0]
                    print(hit.score, best_score)
                    if hit.score > best_score:
                        seen[hit_name] = hit.score, domain, hit

    if seen:
        recombinases = []
        with open(
            outpath / f"{genome_id}.recombinase_scan.tsv",
            "wt",
            encoding="UTF-8",
        ) as rscan_out:
            print(*RECOMBINASE_SCAN_HEADER, sep="\t", file=rscan_out)

            for protein_id, (score, domain, hit) in sorted(seen.items()):
                hmm_name = domain.alignment.hmm_name
                print(protein_id, score, hmm_name)

                recombinase = hmm_name.lower()
                for name, alias in MGE_ALIASES.items():
                    recombinase = recombinase.replace(name, alias)

                rule = mge_rules.get(recombinase)
                if not rule:
                    raise ValueError(f"Cannot find rule for {recombinase}.")

                mges = rule.get_signals()
                confidence = ("low", "high")[len(mges) == 1]

                print(
                    protein_id,
                    recombinase,
                    domain.alignment.hmm_accession,
                    ";".join(mges),
                    hit.evalue,
                    hit.score,
                    confidence,
                    sep="\t",
                    file=rscan_out,
                )

                protein = proteins.get(protein_id)
                if protein is not None:
                    mge_attribs = ";".join(
                        f"{k}={str(v).replace(';', ',')}"
                        for k, v in zip(
                            ("ID", "recombinase", "PFAM", "predicted_mge", "evalue", "score", "confidence",),
                            (protein_id, recombinase, hmm_name, ",".join(mges), hit.evalue, hit.score, confidence,)
                        )
                    )
                    recombinases.append(
                        (
                            protein_id[:protein_id.rfind("_")],
                            "mgexpose",
                            "gene",
                            protein.start,
                            protein.end,
                            f"{hit.score:.5f}",
                            protein.strand,
                            ".",
                            mge_attribs,
                        )
                    )
        print(*recombinases, sep="\n")

        with open(
            outpath / f"{genome_id}.recombinase_scan.gff3",
            "wt",
            encoding="UTF-8",
        ) as rscan_gff:
            print("##gff-version 3", file=rscan_gff)
            for line in sorted(recombinases, key=lambda x: (x[0], int(x[3]), int(x[4]))):
                # gnl|AGSH|NT12270_27_3   dde_tnp_is1     PF03400.12      is_tn   3.1e-74 245.6   high
                print(*line, sep="\t", file=rscan_gff,)
