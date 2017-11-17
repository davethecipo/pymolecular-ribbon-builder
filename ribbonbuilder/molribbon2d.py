#!/usr/bin/python
# -*- coding: utf-8 -*-
import logging
import sys

import openbabel

import ribbonbuilder.translate as t
import ribbonbuilder.obabelwrapper as obwrapper
from ribbonbuilder.parsing import CLIParser, CLIAdapter

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
logger = logging.getLogger(__name__)


def main(args=sys.argv[1:]):
    cli = CLIParser()
    cli.applyCLI(args=args)
    input_file = cli.conf['in']
    oligomer = obwrapper.readFile(input_file)
    polymer = openbabel.OBMol()
    adapter = CLIAdapter(oligomer, cli.conf)
    conf = adapter.conf
    builder = t.MolBuilder(oligomer,
                           polymer,
                           a1=conf['a1'],
                           a2=conf['a2'],
                           n=conf['n'],
                           m=conf['m'],
                           base=conf['base'])
    builder.top_border = conf['top_border']
    builder.bottom_border = conf['bottom_border']
    builder.left_border = conf['left_border']
    builder.right_border = conf['right_border']
    builder.create_all()

    out_file = cli.conf['out']
    obwrapper.writeFile(out_file, polymer)


if __name__ == '__main__':
    main(sys.argv[1:])
