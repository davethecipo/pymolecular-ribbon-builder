#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import openbabel

from ribbonbuilder import translate as t

import pytest

@pytest.fixture()
def base_and_param():
    mol = openbabel.OBMol()
    distance = 3.0
    carbon = 6
    for vec in [(0.0, 0.0, 0.0),
                (distance, 0.0, 0.0),
                (0.0, distance, 0.0),
                (distance, distance, 0.0)]:
        atom = openbabel.OBAtom()
        atom.SetVector(*vec)
        atom.SetAtomicNum(carbon)
        mol.AddAtom(atom)
    return mol, distance

@pytest.fixture()
def empty_mol():
    return openbabel.OBMol()

def test_build_bulk(base_and_param, empty_mol):
    olig, distance = base_and_param
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
