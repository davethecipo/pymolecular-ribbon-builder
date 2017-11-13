#!/usr/bin/python
# -*- coding: utf-8 -*-


from collections import namedtuple
import itertools


class Atom(namedtuple('Atom', 'element index')):
    pass


def _iter_atoms(atom_indices, atom_coordinates, translation, out_atoms):
    """Helper function for build_bulk and build_closure functions"""
    for element, indices in atom_indices.items():
        for index in indices:
            atom_initial_position = atom_coordinates[element][index]
            final_position = atom_initial_position + translation
            # in python 2 dict.items() doesn't reflect changes after the
            # .items() call, but since the dict doesn't change, in this case
            # it doesn't matter
            if element not in out_atoms.keys():
                out_atoms[element] = []
            out_atoms[element].append(final_position)


def build_closure_lower_left(closure_atoms, a1, a2, n, m, atoms, offset):
    """Returns dict of atoms that close the molecule along n2 border.

    We call the lower left border the one where a1=0, a2 ranges from 1 to m.
    closure_atoms: dict containing the atom indices for closing the border,
     for example {'C':[0, 1]} indicates two carbon atoms corresponding to the
     indices 0 and 1
    a1, a2: base vectors
    n, m: integer multiples for repeating a1 and a2 respectively
    atoms: dict containing the atoms positions, like {'C':[np.array([0, 0])]}
     means one carbon atom in the origin
    offset: vector that represents the offset of the closure atoms w.r.t
     the base. For example, once you have a1 and a2 vectors, you can say
     >>> offset = a1*2 + a2*0
     this indicates that atoms are distant 2*a1 from the base"""
    closure_lower_left = {}
    for step_m in range(m):
        tot_translation = step_m*a2 - offset
        _iter_atoms(closure_atoms, atoms, tot_translation, closure_lower_left)
    return closure_lower_left


def build_closure_lower_right(closure_atoms, a1, a2, n, m, atoms, offset):
    closure_lower_right = {}
    for step_n in range(n):
        # since the range function goes from 0 to n-1, I need m-1
        tot_translation = step_n*a1 + (m-1)*a2 - offset
        _iter_atoms(closure_atoms, atoms, tot_translation, closure_lower_right)
    return closure_lower_right


def build_closure_upper_left(closure_atoms, a1, a2, n, m, atoms, offset):
    closure_upper_left = {}
    for step in range(n):
        tot_translation = step*a1 - offset
        _iter_atoms(closure_atoms, atoms, tot_translation, closure_upper_left)
    return closure_upper_left


def build_closure_upper_right(closure_atoms, a1, a2, n, m, atoms, offset):
    closure_upper_right = {}
    for step in range(m):
        tot_translation = (n-1)*a1 + step*a2 - offset
        _iter_atoms(closure_atoms, atoms, tot_translation, closure_upper_right)
    return closure_upper_right


def build_bulk(base_atoms, a1, a2, n, m, atoms):
    """Return a dictionary of the bulk atoms (with coordinates).

    Output example: {'C': [np.array([0, 0]), np.array([1, 1])]}
    In this case there're two carbon atoms in position (0, 0) and (1, 1)
    base_atoms: dict as returned by _build_atoms_dict
    a1 and a2 are the translation base vectors (type numpy.array)
    n is the repetition for a1 and m the one for a2
    atoms is the dict returned by obtain_atoms (atoms with position vectors)
    """
    bulk = {}
    for step_a, step_b in itertools.product(range(n), range(m)):
        tot_translation = step_a*a1 + step_b*a2
        _iter_atoms(base_atoms, atoms, tot_translation, bulk)
    return bulk
        


    
