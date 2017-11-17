#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import openbabel

from ribbonbuilder import translate as t

import pytest


@pytest.fixture()
def mol_and_params():
    mol = openbabel.OBMol()
    distance = 3.0
    border_distance = 1.0
    carbon = 6
    hydrogen = 1
    # Add the 4 "bulk" atoms first
    for vec in [(0.0, 0.0, 0.0),
                (distance, 0.0, 0.0),
                (0.0, distance, 0.0),
                (distance, distance, 0.0)]:
        atom = openbabel.OBAtom()
        atom.SetVector(*vec)
        atom.SetAtomicNum(carbon)
        mol.AddAtom(atom)
    # Now add the atoms on the border
    for vec in [(0.0, -border_distance, 0.0), #  atom 5
                (distance, -border_distance, 0.0), #  atom 6
                (distance+border_distance, 0.0, 0.0), #  atom 7
                (distance+border_distance, distance, 0.0), #  atom 8
                (distance, distance+border_distance, 0.0), #  atom 9
                (0.0, distance+border_distance, 0.0), #  atom 10
                (-border_distance, distance, 0.0), #  atom 11
                (-border_distance, 0.0, 0.0) #  atom 12
                ]:
        atom = openbabel.OBAtom()
        atom.SetVector(*vec)
        atom.SetAtomicNum(hydrogen)
        mol.AddAtom(atom)
    return mol, distance, border_distance


@pytest.fixture()
def border_top_and_param():
    border_dist = 1.0
    atoms = []


@pytest.fixture()
def empty_mol():
    return openbabel.OBMol()


def test_build_bulk(mol_and_params, empty_mol):
    olig, distance, border_distance = mol_and_params
    a1 = np.array([0.0, distance, 0.0])
    a2 = np.array([distance, 0.0, 0.0])
    n = 2
    m = 3
    builder = t.MolBuilder(olig, empty_mol, a1, a2, n, m, [1])
    builder.create_bulk()
    iterator = openbabel.OBMolAtomIter(empty_mol)
    expected = np.array([np.array([0.0, 0.0, 0.0]),
                np.array([distance, 0.0, 0.0]),
                np.array([0.0, distance, 0.0]),
                np.array([distance, distance, 0.0]),
                np.array([0.0, 2*distance, 0.0]),
                np.array([distance, 2*distance, 0.0])])
    result = []
    for atom in iterator:
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        result.append(np.array([x, y, z]))
    result = np.array(result)
    assert len(result) == n*m
    for elem in expected:
        assert elem in result


def test_top(mol_and_params, empty_mol):
    olig, distance, border_distance = mol_and_params
    a1 = np.array([0.0, distance, 0.0])
    a2 = np.array([distance, 0.0, 0.0])
    n = 2
    m = 3
    builder = t.MolBuilder(olig, empty_mol, a1, a2, n, m, [1])
    builder.top_border = {'atoms': [10], 'offset': np.array([0, distance, 0])}
    builder.create_top_border()
    iterator = openbabel.OBMolAtomIter(empty_mol)
    expected = np.array([np.array([0.0, 2*distance+border_distance, 0.0]),
                         np.array([distance, 2*distance+border_distance, 0.0])])
    result = []
    for atom in iterator:
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        result.append(np.array([x, y, z]))
    result = np.array(result)
    assert len(result) == 2
    for elem in expected:
        assert elem in result


def test_bottom(mol_and_params, empty_mol):
    olig, distance, border_distance = mol_and_params
    a1 = np.array([0.0, distance, 0.0])
    a2 = np.array([distance, 0.0, 0.0])
    n = 2
    m = 3
    builder = t.MolBuilder(olig, empty_mol, a1, a2, n, m, [1])
    builder.bottom_border = {'atoms': [5], 'offset': np.array([0, 0, 0])}
    builder.create_bottom_border()
    iterator = openbabel.OBMolAtomIter(empty_mol)
    expected = np.array([np.array([0, -border_distance, 0.0]),
                         np.array([distance, -border_distance, 0.0])])
    result = []
    for atom in iterator:
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        result.append(np.array([x, y, z]))
    result = np.array(result)
    assert len(result) == 2
    for elem in expected:
        assert elem in result


def test_left(mol_and_params, empty_mol):
    olig, distance, border_distance = mol_and_params
    a1 = np.array([0.0, distance, 0.0])
    a2 = np.array([distance, 0.0, 0.0])
    n = 2
    m = 3
    builder = t.MolBuilder(olig, empty_mol, a1, a2, n, m, [1])
    builder.left_border = {'atoms': [11],
                             'offset': np.array([0, distance, 0])}
    builder.create_left_border()
    iterator = openbabel.OBMolAtomIter(empty_mol)
    expected = np.array([np.array([-border_distance, 0, 0]),
                         np.array([-border_distance, distance, 0]),
                         np.array([-border_distance, 2*distance, 0])
                         ])
    result = []
    for atom in iterator:
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        result.append(np.array([x, y, z]))
    result = np.array(result)
    assert len(result) == 3
    for elem in expected:
        assert elem in result


def test_right(mol_and_params, empty_mol):
    olig, distance, border_distance = mol_and_params
    a1 = np.array([0.0, distance, 0.0])
    a2 = np.array([distance, 0.0, 0.0])
    n = 2
    m = 3
    builder = t.MolBuilder(olig, empty_mol, a1, a2, n, m, [1])
    builder.right_border = {'atoms': [8],
                             'offset': np.array([distance, distance, 0])}
    builder.create_right_border()
    iterator = openbabel.OBMolAtomIter(empty_mol)
    expected = np.array([np.array([distance+border_distance, 0, 0]),
                         np.array([distance+border_distance, distance, 0]),
                         np.array([distance+border_distance, 2*distance, 0])
                         ])
    result = []
    for atom in iterator:
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        result.append(np.array([x, y, z]))
    result = np.array(result)
    assert len(result) == 3
    for elem in expected:
        assert elem in result


def test_build_all(mol_and_params, empty_mol):
    olig, distance, border_distance = mol_and_params
    a1 = np.array([0.0, distance, 0.0])
    a2 = np.array([distance, 0.0, 0.0])
    n = 2
    m = 3
    builder = t.MolBuilder(olig, empty_mol, a1, a2, n, m, [1])
    builder.top_border = {'atoms': [10], 'offset': np.array([0, distance, 0])}
    builder.bottom_border = {'atoms': [5], 'offset': np.array([0, 0, 0])}
    builder.left_border = {'atoms': [11],
                           'offset': np.array([0, distance, 0])}
    builder.right_border = {'atoms': [8],
                            'offset': np.array([distance, distance, 0])}
    builder.create_all()
    iterator = openbabel.OBMolAtomIter(empty_mol)
    expected = np.array([np.array([0.0, 0.0, 0.0]), # 1; start bulk
                np.array([distance, 0.0, 0.0]), # 2
                np.array([0.0, distance, 0.0]), # 3
                np.array([distance, distance, 0.0]), # 4
                np.array([0.0, 2*distance, 0.0]), # 5
                np.array([distance, 2*distance, 0.0]), # 6; end of bulk
                np.array([0, -border_distance, 0.0]), # 7; start bottom
                np.array([distance, -border_distance, 0.0]), # 8; end bottom
                np.array([distance + border_distance, 0, 0]), # 9; start right
                np.array([distance + border_distance, distance, 0]), # 10
                np.array([distance + border_distance, 2 * distance, 0]), # 11; end right
                np.array([distance, 2 * distance + border_distance, 0.0]), # 12; start top
                np.array([0.0, 2 * distance + border_distance, 0.0]), # 13; end top
                np.array([-border_distance, 2 * distance, 0]),  # 14; start left
                np.array([-border_distance, distance, 0]),  # 15
                np.array([-border_distance, 0, 0])  # 16; end left
                        ])
    result = []
    for atom in iterator:
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        result.append(np.array([x, y, z]))
    result = np.array(result)
    assert len(result) == 16
    for elem in expected:
        assert elem in result