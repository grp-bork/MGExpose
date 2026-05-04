# pylint: disable=R0903

""" Module contains various reader/parser functions """

import csv
import json
import re
import sys

from ..chunk_reader import get_lines_from_chunks
from ...genes.gene import Gene
from ...islands.genomic_island import GenomicIsland
from ...recombinases import MgeRule
from ...rules.secretion import CONJSCAN_RULES


def read_preannotated_genes(f):
    """ Read genes from previous run via gene_info.txt.
    Returns Gene objects via generator.
    """
    header = None
    for line in get_lines_from_chunks(f):
        line = line.strip().split("\t")
        if header is None:
            header = line
        else:
            line = [(item, None)[item == "None"] for item in line]
            yield Gene.from_geneinfo(**dict(zip(header, line)))


def read_fasta(f):
    """ Read a fasta file. """
    header, seq = None, []
    for line in get_lines_from_chunks(f):
        if line[0] == ">":
            if seq:
                yield header, "".join(seq)
                seq.clear()
            header = line.strip()[1:]
        else:
            seq.append(line.strip())
    if seq:
        yield header, "".join(seq)


def read_recombinase_hits(f):
    """ Read hmmer output from recombinase scan.

    Returns (gene_id, mge_name) tuples via generator.
    """
    
    best_hits = {}

    with open(f, "rt", encoding="UTF-8") as _in:
        lines = _in.readlines()

        fmt = None
        if lines[0].startswith("##gff-version 3"):
            fmt = "gff3"
        elif lines[0].startswith("#unigene"):
            fmt = "mge_pred"
        elif lines[0][0] == "#" and line[0].endswith("domain number estimation ----"):
            fmt = "raw"
        elif lines[0][0] != "#":
            fmt = "best"

        for line in lines:
            line = line.strip()
            if line and line != "#":
                if fmt == "gff3":
                    attribs = dict(item.split("=") for item in line.split("\t")[8].split(";"))
                    gene_id, mge = attribs.get("ID"), attribs.get("recombinase")
                    best_hits[gene_id] = (0.0, mge)
                elif fmt == "mge_pred":
                    gene_id, mge = line.split("\t")[:2]
                    best_hits[gene_id] = (0.0, mge)
                elif fmt == "raw":
                    gene_id, _, mge, _, _, score, *_ = re.split(r"\s+", line)
                    score = float(score)
                    seen_score = best_hits.setdefault(gene_id, [0.0, ""])[0]
                    if score > seen_score:
                        best_hits[gene_id] = (score, mge)
                else:
                    gene_id, mge = line.split("\t")[:2]
                    best_hits[gene_id] = (0.0, mge)
        
        for gene_id, (_, mge) in best_hits.items():
            yield gene_id, mge


        # for line in _in:
        #     line = line.strip()
        #     if line and line[0] != "#":
        #         if pyhmmer:
        #             gene_id, mge = line.split("\t")[:2]
        #         else:
        #             gene_id, _, mge, *_ = re.split(r"\s+", line)
        #         yield gene_id, mge


# would love to add raw scan parsing to annotator,
# but then the upstream filtering doesn't work anymore... >:(
# def read_recombinase_scan(f):
# 	recombinase_hits = {}
# 	with open(f, "rt") as _in:
# 		for line in _in:
# 			line = line.strip()
# 			if line and line[0] != "#":
# 				gene_id, _, mge, pfam_acc, evalue, score, *_ = re.split(r"\s+", line)
# 				score = float(score)
# 				best_hit = recombinase_hits.get(gene_id)
# 				if best_hit is None or score > best_hit[0]:
# 					recombinase_hits[gene_id] = score, mge, pfam_acc, evalue

# 	for gene_id, recombinase_annotation in recombinase_hits.items():
# 		yield gene_id, recombinase_annotation


# def parse_macsyfinder_rules(f, macsy_version=2):
#     """ Read macsyfinder rules.

#     Returns dictionary {secretion_system: {mandatory: count, accessory: count}}.
#     """
#     key_col, mandatory_col, accessory_col = (0, 1, 2) if macsy_version == 2 else (1, 5, 6)

#     with open(f, "rt", encoding="UTF-8") as _in:
#         return {
#             row[key_col].replace("_putative", ""): {
#                 "mandatory": int(row[mandatory_col]),
#                 "accessory": int(row[accessory_col]),
#             }
#             for row_index, row in enumerate(csv.reader(_in, delimiter="\t"))
#             if row_index and row and not row[0].startswith("#")
#         }
def parse_macsyfinder_rules(f):
    """ Read macsyfinder rules. """
    with open(f, "rb") as _in:
        return json.load(_in)


def parse_macsyfinder_report(f, f_rules):
    """ Read macsyfinder/conjscan results.

    Returns (gene_id, conjscan_results) tuples via generator.
    """

    if f_rules:
        rules = parse_macsyfinder_rules(f_rules)
    else:        
        rules = CONJSCAN_RULES

    with open(f, "rt", encoding="UTF-8") as _in:
        d = {}
        for line in _in:
            line = line.strip()
            if line and line[0] != "#" and line[:8] != "replicon":
                # replicon is the start of header line
                cols = re.split(r"\s+", line.strip())
                _, hit_id, gene_name, _, model_fqn, _, _, hit_status, *_ = cols
                system = model_fqn.replace("CONJ/", "")
                if system:
                    rule = rules.get(system)
                    if rule is None:
                        print(
                            "WARNING: cannot find txsscan-rule for system:",
                            f"`{system}`",
                            file=sys.stderr,
                        )
                    d.setdefault(hit_id, []).append((gene_name, system, rule, hit_status))

        yield from d.items()


def read_mge_rules(f, for_recombinase_scan=False,):
    """ Read MGE rules.

    Returns dictionary {mge: MgeRule}.
    """
    with open(f, "rt", encoding="UTF-8") as _in:
        rules = {
            row[0].lower(): MgeRule(row[0], *(tuple(map(int, row[1:]))), for_recombinase_scan)
            for i, row in enumerate(csv.reader(_in, delimiter="\t"))
            if i != 0
        }

    return rules


def read_precomputed_islands(genome_id=None, single_island=None, island_file=None,):
    """ Helper function to deal with precomputed regions/islands. """
    precomputed_islands = None
    if single_island and island_file:
        raise ValueError("Both --single_island and --precomputed_islands set.")
    if single_island and not island_file:
        precomputed_islands = [
            GenomicIsland.from_region_string(single_island, genome_id=genome_id,)
        ]
    elif not single_island and island_file:
        with open(island_file, "rt", encoding="UTF-8",) as _in:
            precomputed_islands = [
                GenomicIsland.from_region_string(line, genome_id=genome_id,) for line in _in
            ]

    if precomputed_islands is not None:
        precomputed_islands_by_contig = {}
        for island in precomputed_islands:
            precomputed_islands_by_contig.setdefault(island.contig, []).append(island)
        precomputed_islands = precomputed_islands_by_contig

    return precomputed_islands
