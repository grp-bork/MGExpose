import pathlib

import pyhmmer

from .gffio import read_prodigal_gff
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


def get_protein_coords(gff):
    proteins = {}
    for gene in read_prodigal_gff(gff):
        gene.id = f'{gene.contig}_{gene.id.split("_")[-1]}'
        proteins[gene.id] = gene
    return proteins


def run_pyhmmer(args):
    proteins = get_protein_coords(args.gff)

    if args.mge_rules and pathlib.Path(args.mge_rules).is_file():
        mge_rules = read_mge_rules(args.mge_rules, recombinase_scan=True)
    else:
        raise ValueError("Cannot read mge_rules.")

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
    # filtered_table_out = open(
    #     outpath / f"{args.genome_id}.recombinase_hmmsearch.besthits.out",
    #     "wb"
    # )

    with raw_table_out:  # filtered_table_out:
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
            # hits.write(filtered_table_out, header=write_header)

    if seen:
        recombinases = []
        with open(
            outpath / f"{args.genome_id}.recombinase_scan.tsv",
            "wt",
            encoding="UTF-8",
        ) as rscan_out:
            print(*RECOMBINASE_SCAN_HEADER, sep="\t", file=rscan_out)

            for protein_id, (score, domain, hit) in sorted(seen.items()):
                hmm_name = domain.alignment.hmm_name.decode()
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
                    domain.alignment.hmm_accession.decode(),
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
                            ("recombinase", "PFAM", "predicted_mge", "evalue", "score", "confidence",),
                            (recombinase, hmm_name, ",".join(mges), hit.evalue, hit.score, confidence,)
                        )
                    )
                    # attrib_str = ";".join(f"{item[0]}={item[1]}" for item in protein.attribs.items() if item[1])
                    recombinases.append(
                        (
                            protein_id[:protein_id.rfind("_")],
                            "proMGE_recombinase_scan",
                            "gene",
                            protein.start,
                            protein.end,
                            f"{hit.score:.5f}",
                            protein.strand,
                            ".",
                            # ";".join((mge_attribs, attrib_str,))
                            mge_attribs,
                        )
                    )

        with open(
            outpath / f"{args.genome_id}.recombinase_scan.gff3",
            "wt",
            encoding="UTF-8",
        ) as rscan_gff:
            print("##gff-version 3", file=rscan_gff)
            for line in sorted(recombinases, key=lambda x: (x[0], int(x[3]), int(x[4]))):
                # gnl|AGSH|NT12270_27_3   dde_tnp_is1     PF03400.12      is_tn   3.1e-74 245.6   high
                print(*line, sep="\t", file=rscan_gff,)
