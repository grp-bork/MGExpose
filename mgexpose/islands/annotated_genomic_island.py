# pylint: disable=C0116,C0301,R0902,R0916,R0913,R0917
"""
Data Structures Module

This module is designed to simplify the handling of different genomic sequences,
including but not limited to:

- Genomic Island
- Annotated Genomic Island
- MGE Genomic Island
- Gene

The end product of the pipeline is MGE Genomic Island, an Annotated Genomic Island
consisting of Genes.
It can be saved in a tsv or gff3 format together with its attributes and gene annotations.
The MGE type of each MGE Genomic Island is defined by applying MGE Rule.
"""
import logging

from collections import Counter
from dataclasses import dataclass

from .genomic_island import GenomicIsland


logger = logging.getLogger(__name__)


@dataclass
class AnnotatedGenomicIsland(GenomicIsland):
    '''The following class extends generic Genomic Island with MGE machinery annotations.'''

    phage_count: int = 0
    conj_count: int = 0
    conj_man_count: int = 0
    valid_entire: bool = False
    valid_mandatory: bool = False
    valid_accessory: bool = False

    def __post_init__(self):
        """ Apply annotations. """
        conjugation_systems = {}
        cm_counts = Counter()

        for gene in self.genes:
            self.phage_count += gene.phage is not None
            if gene.conjugation_systems:
                self.conj_count += 1

                has_mandatory_system = False
                for system, rule in zip(gene.conjugation_systems, gene.conjugation_rules):
                    try:
                        _, system = system.split(":")
                    except ValueError:
                        continue
                    if system.split("/")[1].split("_")[0] in ("dCONJ", "T4SS", "MOB",):
                        has_mandatory_system = True
                    if rule is not None:
                        cm_counts[(system, False)] += 1
                        cm_counts[(system, True)] += 1
                        conjugation_systems[system] = rule

                self.conj_man_count += has_mandatory_system

                logger.info("conjugation system: Gene %s -> conj_count = %s", gene.id, self.conj_man_count)
                # self.conj_man_count += (
                #     # gene.conjugation_system.upper()[:4] in ("CONJ", "T4SS",)
                #     gene.conjugation_system.split("_")[0] in ("dCONJ", "T4SS", "MOB",)
                # )
                # if gene.conjugation_rule is not None:
                #     cm_counts[(gene.conjugation_system, False)] += 1
                #     cm_counts[(gene.conjugation_system, True)] += 1

                #     conjugation_systems[gene.conjugation_system] = gene.conjugation_rule

        # TODO: validate if still works
        for system, rule in conjugation_systems.items():
            self.valid_mandatory |= (cm_counts[(system, True)] >= rule["mandatory"] / 2)
            self.valid_accessory |= (cm_counts[(system, False)] >= rule["accessory"] / 2)
            self.valid_entire |= (
                cm_counts[(system, True)] == rule["mandatory"] and
                cm_counts[(system, False)] == rule["accessory"]
            )

    def __str__(self):
        """ String representation. """
        return "\t".join(
            map(
                str, (
                    self.contig,
                    self.start,
                    self.end,
                    len(self),
                    len(self.genes),
                    sum(self.recombinases.values()),
                    0,  # is still necessary?
                    self.phage_count,
                    self.conj_count,
                    self.conj_man_count,
                    self.valid_mandatory,
                    self.valid_entire,
                    ",".join(self.recombinases),
                )
            )
        )
