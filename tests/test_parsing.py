#!/usr/bin/python
# -*- coding: utf-8 -*-
import pytest

from ribbonbuilder.parsing import CLIParser, RequiredParameter, EmptyBaseError

@pytest.fixture()
def required_args():
    return ['-b', '0', '-n', '2', '-m', '3',
            '-a1', '0', '1', '-a2', '0', '2']

@pytest.fixture()
def argparse_required():
    return ['in.xyz', 'out.xyz']

def test_cli_overrides_config(tmpdir, required_args, argparse_required):
    cli_parser = CLIParser()
    p = tmpdir.mkdir("sub").join("hello.txt")
    p.write("n: 5")
    args = argparse_required + required_args + ['-c', str(p)]
    cli_parser.applyCLI(args)
    assert cli_parser.conf['n'] == 2

def test_raise_missing_required_arg(argparse_required):
    cli_parser = CLIParser()
    args = argparse_required
    with pytest.raises(RequiredParameter):
        cli_parser.applyCLI(args)

def test_empty_base(required_args, argparse_required):
    cli_parser = CLIParser()
    required_args.pop(1)
    args = argparse_required + required_args
    with pytest.raises(EmptyBaseError):
        cli_parser.applyCLI(args)