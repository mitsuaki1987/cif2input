#!/usr/bin/python3
import sys
import pymatgen


args = sys.argv

structure = pymatgen.Structure.from_file(args[1])
print(structure)
print(structure.distance_matrix)
print(structure.sites[0])
print(len(structure.get_neighbors(structure.sites[0], 2.6)))
