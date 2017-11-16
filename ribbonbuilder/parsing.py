import argparse
import logging
import os
import sys

import numpy as np
import yaml

from ribbonbuilder.translate import MolBuilder

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
logger = logging.getLogger(__name__)


class ConfigError(Exception):
    def __init__(self, config, msg=None):
        if msg is None:
            msg = "Error with configuration {}".format(config)
        super(ConfigError, self).__init__(msg)
        self.config = config


class EmptyBaseError(ConfigError):
    def __init__(self, config):
        super(EmptyBaseError, self).__init__(
            config, msg="At least one atom is required for the base.")


class RequiredParameter(ConfigError):
    def __init__(self, config, parameter):
        super(RequiredParameter, self).__init__(
            config, msg="Parameter {} is required.".format(parameter))
        self.parameter = parameter

class Config(object):
    """A container for CLI arguments."""

    def __init__(self):
        # values are set to None because we use this object as the namespace
        # for the cli configuration
        self.base = None
        self.left = None
        self.left_offset = None
        self.right = None
        self.right_offset = None
        self.top = None
        self.top_offset = None
        self.bottom = None
        self.bottom_offset = None
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


class CLIParser(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Build a molecular ribbon starting from'
                        ' a small molecule.')
        self._populate_parser_opts()
        self.cli_config = Config()
        self.conf = self.default_config()

    def _populate_parser_opts(self):
        self.parser.add_argument('-b', '--base', nargs='*',
                                 help='atom indexes for the base')

        self.parser.add_argument('-top', nargs='*', dest="top",
                                 help='atom indexes of the top border')
        self.parser.add_argument('-top-offset', nargs=2, type=int,
                                 dest='top_offset')

        self.parser.add_argument('-bottom', nargs='*', dest="bottom",
                                 help='atom indexes of the bottom border')
        self.parser.add_argument('-bottom-offset', nargs=2, type=int,
                                 dest='bottom_offset')

        self.parser.add_argument('-left', nargs='*', dest="left",
                                 help='atom indexes of the left border')
        self.parser.add_argument('-left-offset', nargs=2, type=int,
                                 dest='left_offset')

        self.parser.add_argument('-right', nargs='*', dest="right",
                                 help='atom indexes of the right border')
        self.parser.add_argument('-right-offset', nargs=2, type=int,
                                 dest='right_offset')

        self.parser.add_argument('-a1', nargs=2,
                                 help='a1 translation vector. The first atom'
                                      ' is the start of the vector, the second'
                                      ' is the end.')
        self.parser.add_argument('-a2', nargs=2,
                                 help='a2 translation vector. '
                                      'Same format as a1')
        self.parser.add_argument('-n', type=int, help='repetition for a1')
        self.parser.add_argument('-m', type=int, help='repetition for a2')
        self.parser.add_argument('-c', '--cell',
                                 help='YAML file with the translation info')
        self.parser.add_argument('in', help='The input filename')
        self.parser.add_argument('out', help='The output filename')
        self.parser.add_argument('--debug', action="store_true",
                                 help='Print debug statements')

    @staticmethod
    def default_config():
        return {'base': None,
                'left': None,
                'left_offset': [0, 0],
                'right': None,
                'right_offset': [0, 0],
                'top': None,
                'top_offset': [0, 0],
                'bottom': None,
                'bottom_offset': [0, 0],
                'a1': None,
                'a2': None,
                'n': None,
                'm': None,
                }

    @staticmethod
    def required_opts():
        return ['base', 'a1', 'a2', 'n', 'm']

    def _check_required_ops(self, config):
        for parameter in self.required_opts():
            if config[parameter] is None:
                raise RequiredParameter(parameter)
        if len(config['base']) == 0:
            raise EmptyBaseError(config)

    def _apply_yaml_config(self, base_config):
        if self.cli_config.cell is not None:
            with open(self.cli_config.cell, 'rt') as f:
                base_config.update(yaml.load(f.read()))

    def _apply_cli_args(self, base_config):
        cli_dict = vars(self.cli_config)
        for parameter, value in cli_dict.items():
            if value is not None:
                base_config[parameter] = value

    def applyCLI(self):
        # use a class because namespace can not be a dict
        self.parser.parse_args(namespace=self.cli_config)
        if self.cli_config.debug:
            logger.setLevel(logging.DEBUG)
        # populate base_config with the values from the config file
        self._apply_yaml_config(self.conf)
        # now if CLI arguments exist, they should override the config file
        self._apply_cli_args(self.conf)
        # after merging all options, check for required options
        self._check_required_ops(self.conf)
        logger.debug('After merging all options and checking required'
                     ' arguments, self.config is {}'.format(self.conf))
        # input and output files are required arguments of the ArgumentParser,
        # therefore I do not need to check their presence in self.conf
        # callers may want to check that the input and output files exist
        self.conf['in'] = os.path.join(os.getcwd(), self.conf['in'])
        self.conf['out'] = os.path.join(os.getcwd(), self.conf['out'])


class CLIAdapter(object):
    """Read a CLI configuration and adapt it for MolBuilder."""
    def __init__(self, molecule, cli_conf):
        self.mol = molecule
        self.cli_conf = cli_conf

    def _vector_from_indexes(self, start, end):
        start_atom = self.mol.GetAtom(start)
        end_atom = self.mol.GetAtom(end)
        start_pos = np.array([start_atom.GetX(),
                              start_atom.GetY(),
                              start_atom.GetZ()])
        end_pos = np.array([end_atom.GetX(),
                            end_atom.GetY(),
                            end_atom.GetZ()])
        return end_pos - start_pos

    def _intlst(self, lst):
        return [int(i) for i in lst]

    @property
    def conf(self):
        config = {'n': self.cli_conf['n'],
                  'm': self.cli_conf['m'],
                  'base': self._intlst(self.cli_conf['base']),
                  'left_border': MolBuilder._default_border(),
                  'right_border': MolBuilder._default_border(),
                  'top_border': MolBuilder._default_border(),
                  'bottom_border': MolBuilder._default_border(),
                  }

        a1 = self._vector_from_indexes(int(self.cli_conf['a1'][0]),
                                       int(self.cli_conf['a1'][1]))
        a2 = self._vector_from_indexes(int(self.cli_conf['a2'][0]),
                                       int(self.cli_conf['a2'][1]))
        config['a1'] = a1
        config['a2'] = a2
        if self.cli_conf['left']:
            a1offset = self.cli_conf['left_offset'][0] * a1
            a2offset = self.cli_conf['left_offset'][1] * a2
            config['left_border'] = {'atoms': self._intlst(self.cli_conf['left']),
                                     'offset': a1offset +  a2offset}
        if self.cli_conf['right']:
            a1offset = self.cli_conf['right_offset'][0] * a1
            a2offset = self.cli_conf['right_offset'][1] * a2
            config['right_border'] = {'atoms': self._intlst(self.cli_conf['right']),
                                     'offset': a1offset +  a2offset}
        if self.cli_conf['top']:
            a1offset = self.cli_conf['top_offset'][0] * a1
            a2offset = self.cli_conf['top_offset'][1] * a2
            config['top_border'] = {'atoms': self._intlst(self.cli_conf['top']),
                                     'offset': a1offset +  a2offset}
        if self.cli_conf['bottom']:
            a1offset = self.cli_conf['bottom_offset'][0] * a1
            a2offset = self.cli_conf['bottom_offset'][1] * a2
            config['bottom_border'] = {'atoms': self._intlst(self.cli_conf['bottom']),
                                     'offset': a1offset +  a2offset}
        return config