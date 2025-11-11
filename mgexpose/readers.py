# pylint: disable=R0903

""" Module contains various reader/parser functions """

import csv
import gzip
import re
import sys

from .chunk_reader import get_lines_from_chunks
from .recombinases import MgeRule


def read_prodigal_gff(f):
    """ Prodigal gff output reader.

    Returns (gene_id, gff_line) tuples via generator.
    """
    for line in get_lines_from_chunks(f):
        line = line.strip()
        if line and line[0] != "#":
            line = line.split("\t")
            _id = [
                item.split("=")[1]
                for item in line[8].split(";")
                if item.startswith("ID")
            ][0]
            # gene_id = f"{line[0]}_{_id.split('_')[1]}"
            # yield gene_id, line
            yield _id, line


def read_recombinase_hits(f):
    """ Read hmmer output from recombinase scan.

    Returns (gene_id, mge_name) tuples via generator.
    """
    with open(f, "rt", encoding="UTF-8") as _in:
        for line in _in:
            line = line.strip()
            if line and line[0] != "#":
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


class EggnogReader:
    """
    Class to read and parse Eggnog annotations.
    Currently, phages are detected
    based on the regex signals in the emapper 'description' field
    """
    EMAPPER_FIELDS = {
        "v1": {"cog_fcat": 11, "description": 12},
        "v2.0.0": {"cog_fcat": 20, "description": 21},
        "v2.0.2": {"cog_fcat": 20, "description": 21},
        "v2.1.0": {"cog_fcat": 9, "description": 10},
        "v2.1.2": {"cog_fcat": 6, "description": 7,
                   "seed_eggNOG_ortholog": 1,
                   "seed_ortholog_evalue": 2,
                   "seed_ortholog_score": 3,
                   "eggnog_ogs": 4,
                   "max_annot_lvl": 5,
                   "goes": 9,
                   "ec": 10,
                   "kegg_ko": 11,
                   "kegg_pathway": 12,
                   "kegg_module": 13,
                   "kegg_reaction": 14,
                   "kegg_rclass": 15,
                   "brite": 16,
                   "cazy": 18,
                   "bigg_reaction": 19,
                   "pfam": 20
                   },
    }

    @staticmethod
    def parse_emapper(f, emapper_version="v2.1.2", phage_annotation=None):
        """ Parses emapper annotations output.
        Returns (gene_id,  phage_signal, eggnog) -> (str, boolean, tuple) tuples via generator.
        """
        def filter_record(key, value, row):
            return value < len(row) and row[value] and row[value] != "-" and key != "description"

        emapper_fields = EggnogReader.EMAPPER_FIELDS.get(emapper_version)
        if emapper_fields is None:
            raise ValueError(f"{emapper_version} is an unknown emapper annotation format.")
        if f.endswith(".gz"):
            emapper_stream = gzip.open(f, "rt")
        else:
            emapper_stream = open(f, "rt", encoding="UTF-8")

        with emapper_stream:
            for row in csv.reader(emapper_stream, delimiter="\t"):
                if row and row[0][0] != "#":
                    gene_id = row[0]
                    # Collect non-empty eggnog attributes
                    eggnog_gene_ann = tuple(
                        (key, row[value])
                        for key, value in emapper_fields.items()
                        if filter_record(key, value, row)
                    )
                    phage_signal = None
                    if phage_annotation is not None:
                        # note: freetext is converted to lower case here,
                        # so REs only have to match against lower!
                        eggnog_freetext = re.sub(
                            r"\s", "_", row[emapper_fields["description"]]
                        ).lower()
                        is_phage = phage_annotation.is_phage(
                            eggnog_freetext, row[emapper_fields["eggnog_ogs"]]
                        )
                        phage_signal = (None, eggnog_freetext)[is_phage]
                    yield gene_id, phage_signal, eggnog_gene_ann

