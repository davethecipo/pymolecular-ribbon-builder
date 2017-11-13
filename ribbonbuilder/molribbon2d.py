#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys

import numpy as np
import os
import ribbonbuilder.translate as t
import yaml
from ribbonbuilder import parsing

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
logger = logging.getLogger(__name__)


class Config(object):
    """Just a container for CLI arguments."""

    def __init__(self):
        # values are just set to None because we use this object as the namespace
        # for the cli configuration
        # default values should be set in argparse for CLI particular purposes
        # (like the debug flag) or in the base_config variable down below
        self.base = None
        self.lower_left = None
        self.ll_offset = None
        self.lower_right = None
        self.lr_offset = None
        self.upper_left = None
        self.ul_offset = None
        self.upper_right = None
        self.ur_offset = None
        self.a1 = None
        self.a2 = None
        self.n = None
        self.m = None
        self.comment = None
        self.debug = None

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        return self.__str__()


def obtain_atoms(lines_list):
    """Returns a dictionary of elements.

    Format: {'element_label': [np.array([x0, y0]), ...]}

    lines_list: list where each element is a string (a line with an atom
      from xyz file)
    Example: {'C': [[0, 0], [1, 1]],
              'H': [[2, 2], [3, 3]] }

    In the example, there're two carbon atoms with position (0, 0) and (1, 1)
    and two hydrogen atoms with position (2, 2) and (3, 3). The list elements
    are actually numpy arrays, not simple lists.
    """
    atoms = {}
    for line in lines_list:
        # ignore the z component
        element, x, y = line.rstrip('\n').strip(' ').split(' ')[:-1]
        if element not in atoms.keys():
            atoms[element] = []
        atoms[element].append(np.array([float(x), float(y)]))

    return atoms


def atoms_dict_from_path(filepath):
    with open(filepath, 'rt') as f:
        raw_lines = f.readlines()[2:]
    return obtain_atoms(raw_lines)


def main():
    parser = argparse.ArgumentParser(
        description='Build a molecular ribbon starting from a small molecule.')
    parser.add_argument('-b', '--base', nargs='*',
                        help='list of atoms for the base')

    parser.add_argument('-ll' '--lower-left', nargs='*', dest="lower_left",
                        help='list of atoms that close the lower left side')
    parser.add_argument('-oll', '--lower-left-offset', nargs=2, type=int,
                        dest='ll_offset')

    parser.add_argument('-lr', '--lower-right', nargs='*', dest="lower_right",
                        help='list of atoms that close the lower right side')
    parser.add_argument('-olr', '--lower-right-offset', nargs=2, type=int,
                        dest='lr_offset')

    parser.add_argument('-ul', '--upper-left', nargs='*', dest="upper_left",
                        help='list of atoms that close the upper left side')
    parser.add_argument('-oul', '--upper-left-offset', nargs=2, type=int,
                        dest='ul_offset')

    parser.add_argument('-ur', '--upper_right', nargs='*', dest='upper_right',
                        help='list of atoms that close the upper right side')
    parser.add_argument('-our', '--upper-right-offset', nargs=2, type=int,
                        dest='ur_offset')

    parser.add_argument('-a1', nargs=2,
                        help='a1 translation vector. The first atom is the start of the vector, the second is the end.')
    parser.add_argument('-a2', nargs=2,
                        help='a2 translation vector. Same format as a1')
    parser.add_argument('-n', type=int, help='repetition for a1')
    parser.add_argument('-m', type=int, help='repetition for a2')
    parser.add_argument('-c', '--cell',
                        help='YAML file containing the translation info')
    parser.add_argument('--comment', help="Comment to put in the output file")
    parser.add_argument('in', help='The input .xyz file')
    parser.add_argument('out', help='The output .xyz I should write to')

    parser.add_argument('--debug', action="store_true",
                        help='Print debug statements')

    cli_config = Config()
    parser.parse_args(namespace=cli_config)  # namespace can't be a dict

    if cli_config.debug:
        logger.setLevel(logging.DEBUG)

    # set here the default values
    base_config = {'base': None,
                   'lower_left': None,
                   'll_offset': [0, 0],
                   'lower_right': None,
                   'lr_offset': [0, 0],
                   'upper_left': None,
                   'ul_offset': [0, 0],
                   'upper_right': None,
                   'ur_offset': [0, 0],
                   'a1': None,
                   'a2': None,
                   'n': None,
                   'm': None,
                   'comment': 'Comment'
                   }
    # populate base_config with the values from the config file
    if cli_config.cell is not None:
        with open(cli_config.cell, 'rt') as f:
            base_config.update(yaml.load(f.read()))  # dictionary
    cli_dict = vars(cli_config)
    # if CLI arguments exist, they should override the config file
    for parameter, value in cli_dict.items():
        if value is not None:
            base_config[parameter] = value
    # Now check for input
    # Note that argparse uses the type hint for n and m, therefore I don't check those
    # the same goes for the offset values, they are type hinted
    # the offset values are also defaulted to 0 in base_config
    # a1 and a2 have the correct length thanks to argparse, don't check those
    required = ['base', 'a1', 'a2', 'n', 'm']
    for parameter in required:
        if base_config[parameter] is None:
            raise Exception(
                'Error: the {} parameter is not set anywhere'.format(
                    parameter))
    if len(base_config['base']) == 0:
        raise Exception("Error: the list of base atoms is empty.")

    logger.debug('base_config is {}'.format(base_config))

    in_file = os.path.join(os.getcwd(), base_config['in'])
    out_file = os.path.join(os.getcwd(), base_config['out'])

    atoms_positions = atoms_dict_from_path(in_file)

    base_list = parsing.parse_cli_atoms(base_config['base'])
    base = parsing.build_atoms_dict(base_list)

    a1_begin, a1_end = parsing.parse_cli_base_vector(base_config['a1'])
    a2_begin, a2_end = parsing.parse_cli_base_vector(base_config['a2'])
    a1 = parsing.base_vector(a1_begin, a1_end, atoms_positions)
    a2 = parsing.base_vector(a2_begin, a2_end, atoms_positions)

    lower_left_list = parsing.parse_cli_atoms(base_config['lower_left'])
    lower_left = parsing.build_atoms_dict(lower_left_list)
    ll_offset = a1 * base_config['ll_offset'][0] + \
                a2 * base_config['ll_offset'][1]

    lower_right_list = parsing.parse_cli_atoms(base_config['lower_right'])
    lower_right = parsing.build_atoms_dict(lower_right_list)
    lr_offset = a1 * base_config['lr_offset'][0] + \
                a2 * base_config['lr_offset'][1]

    upper_left_list = parsing.parse_cli_atoms(base_config['upper_left'])
    upper_left = parsing.build_atoms_dict(upper_left_list)
    ul_offset = a1 * base_config['ul_offset'][0] + \
                a2 * base_config['ul_offset'][1]

    upper_right_list = parsing.parse_cli_atoms(base_config['upper_right'])
    upper_right = parsing.build_atoms_dict(upper_right_list)
    ur_offset = a1 * base_config['ur_offset'][0] + \
                a2 * base_config['ur_offset'][1]

    n = base_config['n']
    m = base_config['m']

    bulk_atoms = t.build_bulk(base, a1, a2, n, m, atoms_positions)
    lower_left_atoms = t.build_closure_lower_left(lower_left, a1, a2, n, m,
                                                  atoms_positions,
                                                  ll_offset)
    lower_right_atoms = t.build_closure_lower_right(lower_right, a1, a2, n, m,
                                                    atoms_positions,
                                                    lr_offset)
    upper_left_atoms = t.build_closure_upper_left(upper_left, a1, a2, n, m,
                                                  atoms_positions,
                                                  ul_offset)
    upper_right_atoms = t.build_closure_upper_right(upper_right, a1, a2, n, m,
                                                    atoms_positions,
                                                    ur_offset)

    with open(out_file, 'wt') as f:
        lines = parsing.create_xyz(base_config['comment'],
                                   bulk_atoms,
                                   lower_left_atoms,
                                   lower_right_atoms,
                                   upper_left_atoms,
                                   upper_right_atoms)
        f.write("\n".join(lines))


if __name__ == '__main__':
    main()
