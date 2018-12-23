#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import logging
import os
import sys

import openbabel

import ribbonbuilder.translate as t
import ribbonbuilder.obabelwrapper as obwrapper
from ribbonbuilder.parsing import CLIParser, CLIAdapter

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
logger = logging.getLogger(__name__)


def defaultMoleculeAbsolutePath():
    thisFileDir =  os.path.dirname(os.path.realpath(__file__))
    testDir = os.path.join(thisFileDir, os.pardir, 'tests')
    testFile = os.path.join(testDir, 'example-input.xyz')
    return os.path.realpath(testFile)


class CLIArgsForTesting(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Minimal CLI arguments for testing purposes throw away code.')
        self._populate_parser_opts()

    def _populate_parser_opts(self):
        self.parser.add_argument('--in',
            help='Input file readable by OpenBabel', 
            default=defaultMoleculeAbsolutePath()
            )

    def arguments(self):
        return vars(self.parser.parse_args())


def moleculeFromCli(args=sys.argv[1:]):
    cli = CLIArgsForTesting()
    args = cli.arguments()
    input_file = args['in']
    molecule = obwrapper.readFile(input_file)
    return molecule
