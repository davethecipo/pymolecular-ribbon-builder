#!/usr/bin/python
# -*- coding: utf-8 -*-
import itertools

import numpy as np
import openbabel


class MolBuilder(object):
    def __init__(self, oligomer, polymer, a1, a2, n, m, base):
        self.oligomer = oligomer
        self.polymer = polymer
        self.a1 = a1
        self.a2 = a2
        self.n = n
        self.m = m
        self.top_border = self._default_border()
        self.bottom_border = self._default_border()
        self.left_border = self._default_border()
        self.right_border = self._default_border()
        self.base = base

    @staticmethod
    def _default_border():
        return {'atoms': [],
                'offset': np.ndarray(shape=(3,1),
                                     dtype=float,
                                     buffer=np.array([0, 0, 0])
                                     )
                }

    def create_bulk(self):
        for step_a, step_b in itertools.product(range(self.n), range(self.m)):
            tot_translation = step_a * self.a1 + step_b * self.a2
            self._iter_atoms(self.base, tot_translation)

    def create_left_border(self):
        for step_m in range(self.m):
            translation = step_m * self.a2 - self.left_border['offset']
            self._iter_atoms(self.left_border['atoms'], translation)

    def create_right_border(self):
        for step_m in range(self.m):
            # since the range function goes from 0 to m-1, I need m-1
            translation = (self.n -1) * self.a1 + step_m * self.a2 - \
                              self.right_border['offset']
            self._iter_atoms(self.right_border['atoms'], translation)

    def create_top_border(self):
        for step_n in range(self.n):
            translation = step_n * self.a1 + (self.m - 1) * self.a2 - \
                          self.top_border['offset']
            self._iter_atoms(self.top_border['atoms'], translation)

    def create_bottom_border(self):
        for step_n in range(self.n):
            translation = step_n * self.a1 - self.bottom_border['offset']
            self._iter_atoms(self.bottom_border['atoms'], translation)

    def create_all(self):
        self.create_bulk()
        self.create_top_border()
        self.create_bottom_border()
        self.create_left_border()
        self.create_right_border()

    def _iter_atoms(self, atom_indexes, translation):
        """Apply translation to all atoms corresponding to atom_indexes.

        Populates self.polymer with the atoms corresponding to atom_indexes
        after applying the translation.
        """
        for index in atom_indexes:
            atom = self.oligomer.GetAtom(index)
            initial_pos = np.array([atom.GetX(), atom.GetY(), atom.GetZ()])
            final_pos = initial_pos + translation
            new_atom = openbabel.OBAtom()
            new_atom.SetAtomicNum(atom.GetAtomicNum())
            new_atom.SetVector(final_pos[0], final_pos[1], final_pos[2])
            self.polymer.AddAtom(new_atom)
