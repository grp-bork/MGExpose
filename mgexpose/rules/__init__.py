""" MGExpose rule module """

import pathlib

from ..utils.readers import read_mge_rules
from .recombinases import MGE_RULES


def get_recombinase_rules(fn, for_recombinase_scan=False,):
    """ Obtain recombinase rules """
    if fn:
        if not pathlib.Path(fn).is_file():
            raise ValueError(f"Cannot read mge_rules from {fn}.")

        return read_mge_rules(fn, for_recombinase_scan=for_recombinase_scan,)


    mge_rules = MGE_RULES
    if not for_recombinase_scan:
        for rule in mge_rules.values():
            rule._tn3_ce_modify()

    return mge_rules
