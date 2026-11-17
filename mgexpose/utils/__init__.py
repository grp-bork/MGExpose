""" module docstring """

import hashlib

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def compute_hash(seq):
    """ computes a hash for an input sequence taking into account both fwd and rc """
    seq = seq.upper()
    seq_rc = seq.translate(COMPLEMENT)[::-1]

    return hashlib.sha256(min(seq, seq_rc).encode("ascii")).hexdigest()
