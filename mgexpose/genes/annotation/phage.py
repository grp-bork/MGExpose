# pylint: disable=R0903

""" Phage detection via keyword search """

import re


class PhageDetection:
    """ Class to detect phage signals in freetext functional gene annotation. """
    # VIRAL_STRUCTURES = re.compile(
    #     "portal|fiber|collar|terminase|prohead"
    #     "|baseplate|sheath|base-plate|tail|head|capsid|tube"
    # )
    # m/portal|tail_fiber|terminase|prohead|baseplate|tail_sheath|Tail_sheath|
    # base-plate|tail_protein|capsid|tail_tube/
    VIRAL_STRUCTURES_KEYWORDS = (
        r"base-?plate",
        r"capsid",
        r"portal",
        r"prohead",
        r"terminase",
        r"tail_(fiber|protein|sheath|tube)",
    )
    VIRAL_STRUCTURES = re.compile(r"|".join(VIRAL_STRUCTURES_KEYWORDS))

    # EXCLUDE_LIST = re.compile(
    #     "ribosome|ribosomal|30s|50s|sipc|tafi"
    #     "|post-translational|mycolic|macrophage"
    # )
    # m/ribosome|ribosomal|30S|50S|SipC|Tafi|post-translational|mycolic|
    # macrophage|peptidoglycan|sickle|Rhophilin|ATPase|myelin|Cysteine/i)
    EXCLUDED_KEYWORDS = (
        r"ribosom(e|al)",
        r"[35]0s",
        r"sipc",
        r"tafi",
        r"post-?translational",
        r"my(elin|colic)",
        r"macrophage",
        r"peptido-?glycan",
        r"sickle",
        r"rhophilin",
        r"atpase",
        r"cysteine",
    )
    EXCLUDE_LIST = re.compile(r"|".join(EXCLUDED_KEYWORDS))

    # EXPECTED_PHAGES = re.compile("phage|bacteriophage|prophage|lamboid|lambda")
    # m/phage|bacteriophage|prophage|lamboid|lambda|\bMu\b|Mu-like
    PHAGE_KEYWORDS = (
        r"(bacterio|pro)?phage",
        r"lamb(da|oid)",
        r"mu(-like)?"
    )
    EXPECTED_PHAGES = re.compile(r"|".join(PHAGE_KEYWORDS))

    # EXTENDED_VIRAL_STRUCTURES = re.compile("holi|dna-packaging|mu-like|lysis|associated|membrane")
    # m/holi|DNA-packaging|portal|fiber|collar|terminase|prohead|baseplate|
    # sheath|base-plate|lysis|membrane_protein|tail|head|capsid|tube
    EXTENDED_VIRAL_STRUCTURES_KEYWORDS = (
        r"holi",
        r"dna-packaging",
        r"portal",
        r"fiber",
        r"terminase",
        r"prohead",
        r"base-?plate",
        r"sheath",
        r"lysis",
        r"membrane_protein",
        r"tail",
        r"head",
        r"capsid",
        r"tube",
    )
    EXTENDED_VIRAL_STRUCTURES = re.compile(r"|".join(EXTENDED_VIRAL_STRUCTURES_KEYWORDS))

    # EXCLUDE_INTEGRASE = re.compile("[pP]hage[ _]integrase")
    # m/Phage integrase|phage integrase/i
    INTEGRASE = re.compile(r"phage[ _]integrase")

    def __init__(self, phage_filter_file=None):
        self.phage_filter = set()
        if phage_filter_file is not None:
            with open(phage_filter_file, "rt", encoding="UTF-8") as _in:
                self.phage_filter = set(line.strip().split("\t")[0] for line in _in)

    def is_phage(self, eggnog_freetext, eggnog_og):
        """ Filters phage-related parsed eggnog mapper output.
        Returns binary phage signal.
        """
        if eggnog_og in self.phage_filter:
            return False

        viral_structure = PhageDetection.VIRAL_STRUCTURES.search(eggnog_freetext)
        excluded_term = PhageDetection.EXCLUDE_LIST.search(eggnog_freetext)
        if viral_structure and not excluded_term:
            return True

        phage_structure = PhageDetection.EXPECTED_PHAGES.search(eggnog_freetext)
        ext_viral_structure = PhageDetection.EXTENDED_VIRAL_STRUCTURES.search(eggnog_freetext)
        integrase = PhageDetection.INTEGRASE.search(eggnog_freetext)

        if phage_structure and ext_viral_structure and not integrase:
            return True

        return False

        # if all((
        #         PhageDetection.VIRAL_STRUCTURES.search(eggnog_freetext),
        #         not PhageDetection.EXCLUDE_LIST.search(eggnog_freetext),
        # )):
        #     self.phage_annotated.add(gene_id)
        #     return True
        # if all((
        #         gene_id not in self.phage_annotated,
        #         PhageDetection.EXPECTED_PHAGES.search(eggnog_freetext),
        #         any((
        #                 PhageDetection.VIRAL_STRUCTURES.search(eggnog_freetext),
        #                 PhageDetection.EXTENDED_VIRAL_STRUCTURES.search(eggnog_freetext),
        #         )),
        #         not PhageDetection.EXCLUDE_INTEGRASE.search(eggnog_freetext),
        # )):
        #     self.phage_annotated.add(gene_id)
        #     return True
        # return False
