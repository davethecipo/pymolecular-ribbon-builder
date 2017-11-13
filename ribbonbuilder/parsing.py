# convert from the format of the CLI to the format needed by the translating functions
# also function to write xyz file
import re

from .translate import Atom

def parse_cli_atom(atom_string, count_from1 = True):
    atom_elem, atom_index = re.findall(r"[^\W\d_]+|\d+", atom_string)
    atom_index = int(atom_index)
    if count_from1:
        atom_index -= 1
    return Atom(atom_elem, atom_index)

def parse_cli_atoms(atoms_list, **kwargs):
    """Returns a list of Atom namedtuples
    
    atoms_list: a list of strings in the format 'ElementX', for example 'C1', 'H2'
    which stand for first carbon atom and second hydrogen atom"""
    out = []
    for atom in atoms_list:
        out.append(parse_cli_atom(atom, **kwargs))
    return out

def parse_cli_base_vector(atoms_list, **kwargs):
    """Return two Atom namedtuples from the list
    
    Example:
        >>> vector_from_cli = ['C1', 'H2']
        >>> start, end = parse_cli_base_vector(vector_from_cli)
        >>> start
        Atom(element='C', index=0)
        >>> end
        Atom(element='H', index=1)

        The reason why the indices are 0 and 1 (instead of 1 and 2) is that
        Avogadro starts counting from 1, while internally we count from 0.
        """
    start_atom, end_atom = parse_cli_atoms(atoms_list, **kwargs)
    return start_atom, end_atom


def base_vector(start, end, atoms_dict):
    """Return a numpy array in the form [x, y] of the base vector components.
    
    start: an Atom namedtuple
    end: same as start
    atoms_dict: dictionary of atoms
    """
    # I can use the numpy vector subtraction directly
    vector = atoms_dict[end.element][end.index] - atoms_dict[start.element][start.index]
    return vector

def build_atoms_dict(atoms_list):
    """Returns a dict in the form {'elem': [index]}

    atoms_list: a list of Atom namedtuple instances
    
    Example: {'C': [0,1,3], 'H': [1,2]}
    Therefore in this case the returned atoms are three carbon atoms
    (corresponding to the indices 0,1 and 3) and two hydrogen atoms
    (corresponding to the indices 1 and 2)
    """
    atoms = {}
    for atom in atoms_list:
        if atom.element not in atoms.keys():
            atoms[atom.element] = []
        atoms[atom.element].append(atom.index)
    return atoms


def create_xyz(comment, *atom_dicts):
    """Return a list of the lines for the xyz file

    atom_dicts: dict as returned by molribbon2d.obtain_atoms
    comment: this is inserted in the 2nd line

    The function molribbon2d.obtain_atoms returns something like this
    >>> {'C': [np.array([0, 0]), np.array([0, 1])]}
    (two carbon atoms, one in the origin and one in 0, 1)
    numpy arrays are not needed, you can also pass python lists
    >>> {'C': [[0, 0], [0, 1]]}
    Example:
    >>> comment = 'comment'
    >>> base = {'C': [[0, 0], [0, 1]]}
    >>> closure = {'H': [[1, 1], [2, 1]]}
    >>> result = create_xyz(comment, base, closure)
    ['4',
    'comment',
    'C 0.0 0.0 0.0',
    'C 0.0 1.0 0.0',
    'H 1.0 1.0 0.0',
    'H 2.0 1.0 0.0']

    """
    out = []
    for atom_dict in atom_dicts: # base or closures
        for element, coords_list in atom_dict.items():
            for coord_pair in coords_list:
                out.append("{} {} {} 0.0".format(
                    element,
                    float(coord_pair[0]),
                    float(coord_pair[1])))
    length = len(out)
    out.insert(0, str(length))
    out.insert(1, comment)
    return out
        
