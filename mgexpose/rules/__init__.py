import pathlib

from ..utils.readers import read_mge_rules
from .recombinases import MGE_RULES


def get_recombinase_rules(fn, for_recombinase_scan=False,):
    if fn:
        if pathlib.Path(fn).is_file():
            return read_mge_rules(fn, for_recombinase_scan=for_recombinase_scan,)
        else:
            raise ValueError(f"Cannot read mge_rules from {fn}.")
    else:
        mge_rules = MGE_RULES
        if not for_recombinase_scan:
            for rule in mge_rules.values():
                rule._tn3_ce_modify()

        return mge_rules
