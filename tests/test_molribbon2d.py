#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

import openbabel
from scripttest import TestFileEnvironment

from ribbonbuilder import molribbon2d

def test_write_output_cli():
    testdir = os.path.abspath(os.path.dirname(__file__))
    input_fname = os.path.join(testdir, 'example-input.xyz')
    output_fname = os.path.join(testdir, 'example-output.xyz')
    cell_fname = os.path.join(testdir, 'example-cell.yaml')
    base_folder = './scratch'
    env = TestFileEnvironment(base_folder)
    result_fname = 'result.xyz'
    result = env.run('molribbon2d -c {} -n 2 -m 3 {} {}'.format(cell_fname,
                                                                input_fname,
                                                                result_fname))
    assert result_fname in result.files_created


def test_write_correct_output_cli():
    testdir = os.path.abspath(os.path.dirname(__file__))
    input_fname = os.path.join(testdir, 'example-input.xyz')
    output_fname = os.path.join(testdir, 'example-output.xyz')
    cell_fname = os.path.join(testdir, 'example-cell.yaml')
    base_folder = './scratch'
    env = TestFileEnvironment(base_folder)
    result_fname = 'result.xyz'
    result = env.run('molribbon2d -c {} -n 2 -m 3 {} {}'.format(cell_fname,
                                                                input_fname,
                                                                result_fname))
    result_abs_fname = os.path.join(base_folder, result_fname)
    result_mol = openbabel.OBMol()
    expected_mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('xyz', 'xyz')
    conv.ReadFile(expected_mol, output_fname)
    conv.ReadFile(result_mol, result_abs_fname)
    num_atoms = expected_mol.NumAtoms()
    expected_vecs = []
    result_vecs = []
    for idx in range(1, num_atoms+1):
        # TODO this is brittle because it depends on atom indexes which
        #  according to openbabel docs will change. It can be used for now
        expected_atom = expected_mol.GetAtom(idx)
        expected_vec = (expected_atom.GetX(), expected_atom.GetY(),
                        expected_atom.GetZ())
        expected_vecs.append(expected_vec)
        result_atom = result_mol.GetAtom(idx)
        result_vec = (result_atom.GetX(), result_atom.GetY(), result_atom.GetZ())
        result_vecs.append(result_vec)
    for vec in expected_vecs:
        assert vec in result_vecs


def test_write_correct_result_module(tmpdir):
    testdir = os.path.abspath(os.path.dirname(__file__))
    input_fname = os.path.join(testdir, 'example-input.xyz')
    output_fname = os.path.join(testdir, 'example-output.xyz')
    cell_fname = os.path.join(testdir, 'example-cell.yaml')
    tmpfile = str(tmpdir.mkdir('sub').join('result.xyz'))

    n, m = (2, 3)
    molribbon2d.main([input_fname, tmpfile, '-c', cell_fname, '-n', str(n),
                      '-m', str(m)])
    result_mol = openbabel.OBMol()
    expected_mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('xyz', 'xyz')
    conv.ReadFile(expected_mol, output_fname)
    conv.ReadFile(result_mol, tmpfile)
    num_atoms = expected_mol.NumAtoms()
    expected_vecs = []
    result_vecs = []
    for idx in range(1, num_atoms+1):
        # TODO this is brittle because it depends on atom indexes which
        #  according to openbabel docs will change. It can be used for now
        expected_atom = expected_mol.GetAtom(idx)
        expected_vec = (expected_atom.GetX(), expected_atom.GetY(),
                        expected_atom.GetZ())
        expected_vecs.append(expected_vec)
        result_atom = result_mol.GetAtom(idx)
        result_vec = (result_atom.GetX(), result_atom.GetY(), result_atom.GetZ())
        result_vecs.append(result_vec)
    for vec in expected_vecs:
        assert vec in result_vecs