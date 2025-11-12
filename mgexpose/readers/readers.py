# pylint: disable=R0903

""" Module contains various reader/parser functions """

import csv
import gzip
import re
import sys

from ..utils.chunk_reader import get_lines_from_chunks
from ..gene import Gene
from ..recombinases import MgeRule



def read_fasta(f):
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


def read_recombinase_hits(f, pyhmmer=True):
    """ Read hmmer output from recombinase scan.

    Returns (gene_id, mge_name) tuples via generator.
    """
    with open(f, "rt", encoding="UTF-8") as _in:
        for line in _in:
            line = line.strip()
            if line and line[0] != "#":
                if pyhmmer:
                    gene_id, mge = line.split("\t")[:2]
                else:
                    gene_id, _, mge, *_ = re.split(r"\s+", line)
                yield gene_id, mge


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


def parse_macsyfinder_rules(f, macsy_version=2):
    """ Read macsyfinder rules.

    Returns dictionary {secretion_system: {mandatory: count, accessory: count}}.
    """
    key_col, mandatory_col, accessory_col = (0, 1, 2) if macsy_version == 2 else (1, 5, 6)

    with open(f, "rt", encoding="UTF-8") as _in:
        return {
            row[key_col].replace("_putative", ""): {
                "mandatory": int(row[mandatory_col]),
                "accessory": int(row[accessory_col]),
            }
            for row_index, row in enumerate(csv.reader(_in, delimiter="\t"))
            if row_index and row and not row[0].startswith("#")
        }


def parse_macsyfinder_report(f, f_rules, macsy_version=2):
    """ Read macsyfinder/txsscan results.

    Returns (gene_id, txsscan_results) tuples via generator.
    """

    rules = parse_macsyfinder_rules(f_rules, macsy_version=macsy_version)

    key_col, col1, col2 = (1, 4, 8) if macsy_version == 2 else (0, 6, 9)

    with open(f, "rt", encoding="UTF-8") as _in:
        for line in _in:
            line = line.strip()
            if line and line[0] != "#":
                line = re.split(r"\s+", line.strip())
                system = line[col1].replace("TXSS/", "")
                rule = rules.get(system)
                if rule is None:
                    print(
                        "WARNING: cannot find txsscan-rule for system:",
                        f"`{system}`",
                        file=sys.stderr,
                    )

                if line and line[0] and line[0] != "replicon":
                    yield line[key_col], (system, rule, line[col2])


def read_mge_rules(f, recombinase_scan=False):
    """ Read MGE rules.

    Returns dictionary {mge: MgeRule}.
    """
    with open(f, "rt", encoding="UTF-8") as _in:
        rules = {
            row[0].lower(): MgeRule(row[0], *(tuple(map(int, row[1:]))), recombinase_scan)
            for i, row in enumerate(csv.reader(_in, delimiter="\t"))
            if i != 0
        }

    # #special case for Tn3 since it can carry conjugative system#
    # for rule_id, rule in rules.items():
    # 	if "tn3" in rule_id:
    # 		rule.ce = 1

    return rules
