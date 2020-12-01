#!/usr/bin/python3
import sys
import pymatgen
import numpy

args = sys.argv

print("CRYSTAL")
print("CELLP")
structure = pymatgen.Structure.from_file(args[1])
lattice = structure.lattice
print(lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta, lattice.gamma)
print(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

print("STRUC")
icount = 0
for site1 in structure.sites:
    icount += 1
    frac1 = numpy.array([site1.a, site1.b, site1.c])
    print(icount, site1.species_string, site1.species_string, 1.0, frac1[0], frac1[1], frac1[2], "1a", 1)
    print(0.0, 0.0, 0.0, 0.0)
    sites = structure.get_neighbors(site1, 2.6)
print(0, 0, 0, 0, 0, 0, 0)

print("VECTR")
icount = 0
jcount = 0
for site1 in structure.sites:
    icount += 1
    frac1 = numpy.array([site1.a, site1.b, site1.c])
    sites = structure.get_neighbors(site1, 2.6)
    for site2 in sites:
        jcount += 1
        frac2 = numpy.array([site2.a, site2.b, site2.c])
        dfrac = frac2[:] - frac1[:]
        dcart = structure.lattice.get_cartesian_coords(dfrac)
        print(jcount, dcart[0], dcart[1], dcart[2], 0)
        print(icount, 0, 0, 0, 0)
        print(0, 0, 0, 0, 0)
print(0, 0, 0, 0, 0)

print("VECTT")
jcount = 0
for site1 in structure.sites:
    sites = structure.get_neighbors(site1, 2.6)
    for site2 in sites:
        jcount += 1
        print(jcount, 0.5, 255, 0, 0, 0)
print("0 0 0 0 00 0 0 0 0")
