from itertools import combinations

from .obabelwrapper import readFile
from .cliboilerplate import moleculeFromCli
from .maindirections import growthDirectionsOf, biggestDirectionFrom, positionVectorsOf

from matplotlib import pyplot as plt

def main():
    molecule = moleculeFromCli()

    positions = positionVectorsOf(molecule)

    directions = growthDirectionsOf(molecule)
    biggest = biggestDirectionFrom(directions)

    conns = allConnectionsAmong(positions)
    print(conns)


def allConnectionsAmong(positions):
    segments = []
    for combination in combinations(positions, 2):
        segment = combination[1] - combination[0]
        segments.append(segment)
    return segments

if __name__ == "__main__":
    main()