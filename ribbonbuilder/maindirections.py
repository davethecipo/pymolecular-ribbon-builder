import numpy as np
import openbabel

from numpy.linalg import eigh


def mainGrowthDirectionOf(molecule):
    directions = growthDirectionsOf(molecule)
    return biggestDirectionFrom(directions)


def biggestDirectionFrom(directions):
    norms = [np.linalg.norm(i) for i in directions]
    biggestNormIndex = np.argmax(norms)
    return directions[biggestNormIndex]


def growthDirectionsOf(molecule):
    positions = positionVectorsOf(molecule)    
    directions = mainDirectionsFrom(positions)
    return directions


def positionVectorsOf(molecule):
    positions = []

    for atom in openbabel.OBMolAtomIter(molecule):
        position = np.array([atom.GetX(), atom.GetY(), atom.GetZ()])
        positions.append(position)

    return np.array(positions)


def mainDirectionsFrom(positions):
    covarianceMatrix = np.cov(positions, rowvar=False)
    eigenvalues, eigenvectors = eigh(covarianceMatrix)

    directions = []
    for index, eigenvalue in enumerate(eigenvalues):
        direction = eigenvectors[:,index] * eigenvalue
        directions.append(direction)

    return directions
