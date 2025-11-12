import pathlib

import pyhmmer

from .recombinases import MGE_ALIASES
from .readers import read_mge_rules


RECOMBINASE_SCAN_HEADER = (
    "#unigene",
    "recombinase_SMART_hmm_name",
    "PFAM_accession",
    "MGE_prediction",
    "hmmsearch_fullsequence_evalue",
    "hmmsearch_fullsequence_score",
    "MGE_prediction_confidence",
)


def run_pyhmmer(args):
    with pyhmmer.easel.SequenceFile(args.proteins_fasta, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
        protein_seqs = list(seq_file)
    with pyhmmer.plan7.HMMFile(args.recombinase_hmms) as hmm_file:
        hmm_hits = list(
            pyhmmer.hmmsearch(hmm_file, protein_seqs, cpus=args.threads, backend="multiprocessing", bit_cutoffs="gathering")
        )

    outpath = pathlib.Path(args.output_dir)
    outpath.mkdir(exist_ok=True, parents=True,)

    raw_table_out = open(
        outpath / f"{args.genome_id}.recombinase_hmmsearch.out",
        "wb"
    )
    filtered_table_out = open(
        outpath / f"{args.genome_id}.recombinase_hmmsearch.besthits.out",
        "wb"
    )

    with raw_table_out, filtered_table_out:
        seen = {}
        for i, hits in enumerate(hmm_hits):
            write_header = i == 0
            hits.write(raw_table_out, header=write_header)
            for hit in hits:
                hit_name = hit.name.decode()
                for domain in hit.domains:
                    best_score = seen.setdefault(hit_name, (0.0, None, None))[0]
                    print(hit.score, best_score)
                    if hit.score > best_score:
                        seen[hit_name] = hit.score, domain, hit
            hits.write(filtered_table_out, header=write_header)

    if seen and args.mge_rules and pathlib.Path(args.mge_rules).is_file():
        mge_rules = read_mge_rules(args.mge_rules, recombinase_scan=True)

        with open(
            outpath / f"{args.genome_id}.recombinase_scan.tsv",
            "wt",
            encoding="UTF-8",
        ) as rscan_out:
            print(*RECOMBINASE_SCAN_HEADER, sep="\t", file=rscan_out)

            for protein, (score, domain, hit) in seen.items():
                hmm_name = domain.alignment.hmm_name.decode()
                print(protein, score, hmm_name)

                recombinase = hmm_name.lower()
                for name, alias in MGE_ALIASES.items():
                    recombinase = recombinase.replace(name, alias)

                rule = mge_rules.get(recombinase)
                if not rule:
                    raise ValueError(f"Cannot find rule for {recombinase}.")

                mges = rule.get_signals()
                confidence = ("low", "high")[len(mges) == 1]

                print(
                    protein,
                    recombinase,
                    domain.alignment.hmm_accession.decode(),
                    ";".join(mges),
                    hit.evalue,
                    hit.score,
                    confidence,
                    sep="\t",
                    file=rscan_out,
                )