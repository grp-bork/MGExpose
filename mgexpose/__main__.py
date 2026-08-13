#!/usr/bin/env python

# pylint: disable=R0912,R0914,R0915,R0913,R0917

""" Mobile genetic element annotation """

import logging
import sys

from .handle_args import handle_args
from .modules.denovo import denovo
from .modules.gene_calling import gene_calling
from .modules.reannotate import reannotate
from .modules.recombinase_scan import recombinase_scan


logger = logging.getLogger(__name__)


def main():
    """ main """

    args = handle_args(sys.argv[1:])
    logger.info("ARGS: %s", str(args))

    if args.command == "denovo":
        denovo(args)

    elif args.command == "reannotate":
        reannotate(args)

    elif args.command == "call_genes":
        gene_calling(args)

    elif args.command == "recombinase_scan":
        recombinase_scan(args)

    else:
        raise NotImplementedError(f"{args.command} module not implemented.")


if __name__ == "__main__":
    main()
