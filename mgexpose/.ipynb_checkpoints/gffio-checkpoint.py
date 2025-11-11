from .gene import Gene
from .islands import GenomicIsland, MgeGenomicIsland

from .base_logger import logger

def read_genomic_islands_gff(fn):
    with open(fn, "rt", encoding="UTF-8") as _in:
        island = None
        for line in _in:
            line = line.strip()
            if line and line[0] != "#":
                cols = line.split("\t")
                if cols[2] == "region":
                    if island is not None:
                        yield island
                    island = GenomicIsland.from_gff(*cols)
                elif cols[2] == "gene":
                    gene = Gene.from_gff(*cols)
                    if island is not None:
                        island.genes.add(gene)
                    else:
                        raise ValueError("Found gene but no island.")
        if island is not None:
            yield island
       
def read_mge_genomic_islands_gff(fn, relevant_ids=None):
    """
    Generator function to read and parse MGEs and associated genes from a GFF file.

    Parameters:
    - fn: Path to the GFF file.
    - relevant_ids: Optional set of relevant MGE IDs to filter. If None, all MGEs are processed.

    Yields:
    - MgeGenomicIsland objects that match the relevant IDs or all if None.
    """
    with open(fn, "rt", encoding="UTF-8") as _in:
        island = None
        for line in _in:
            line = line.strip()
            if line and line[0] != "#":
                cols = line.split("\t")
                attributes = {kv.split('=')[0]: kv.split('=')[1] for kv in cols[8].split(';') if '=' in kv}

                if cols[2] == "mobile_genetic_element":
                    mge_id = attributes.get("ID")

                    if relevant_ids is None or mge_id in relevant_ids:
                        if island is not None:
                            yield island
                        island = MgeGenomicIsland.from_gff(*cols)

                elif cols[2] == "gene":
                    parent_id = attributes.get("Parent")

                    if island is not None:
                        if relevant_ids is None or parent_id in relevant_ids:
                            gene = Gene.from_gff(*cols)
                            island.genes.add(gene)
                        else:
                            continue
                    else:
                        # This situation should not happen unless the GFF is malformed
                        raise ValueError("Found gene with no preceding island.")
                        
        if island is not None:
            yield island
