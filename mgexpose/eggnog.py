import csv
import gzip
import re


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
