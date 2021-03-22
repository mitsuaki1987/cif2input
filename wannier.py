#!/usr/bin/python3
import sys
import pymatgen
import numpy
import math

args = sys.argv

ntheta = 30
nbase = ntheta**3

base = numpy.ndarray((nbase, 3, 3), dtype=float)
basex = numpy.ndarray((3, 3), dtype=float)
basey = numpy.ndarray((3, 3), dtype=float)
basez = numpy.ndarray((3, 3), dtype=float)

ibase = 0
for ix in range(ntheta):
    thetax = 2.0 * math.pi * ix / ntheta
    basex[0:3, 0:3] = [[1.0,              0.0,               0.0],
                       [0.0, math.cos(thetax), -math.sin(thetax)],
                       [0.0, math.sin(thetax),  math.cos(thetax)]]
    for iy in range(ntheta):
        thetay = 2.0 * math.pi * iy / ntheta
        basey[0:3, 0:3] = [[math.cos(thetay),  0.0, math.sin(thetay)],
                           [0.0,               1.0,              0.0],
                           [-math.sin(thetay), 0.0, math.cos(thetay)]]
        for iz in range(ntheta):
            thetaz = 2.0 * math.pi * iz / ntheta
            basez[0:3, 0:3] = [[math.sin(thetaz),  math.cos(thetaz), 0.0],
                               [math.cos(thetaz), -math.sin(thetaz), 0.0],
                               [0.0,                            0.0, 1.0]]
            base[ibase, 0:3, 0:3] = numpy.dot(numpy.dot(basex, basey), basez)
            ibase += 1

with open(args[1][:-3] + "vesta", 'w') as f:

    print("CRYSTAL", file=f)
    print("CELLP", file=f)
    structure = pymatgen.Structure.from_file(args[1])
    lattice = structure.lattice
    print(lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta, lattice.gamma, file=f)
    print(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, file=f)

    print("STRUC", file=f)
    icount = 0
    for site1 in structure.sites:
        icount += 1
        frac1 = numpy.array([site1.a, site1.b, site1.c])
        print(icount, site1.species_string, site1.species_string, 1.0,
              frac1[0], frac1[1], frac1[2], "1a", 1, file=f)
        print(0.0, 0.0, 0.0, 0.0, file=f)
    print(0, 0, 0, 0, 0, 0, 0, file=f)

    print("VECTR", file=f)
    icount = 0
    for site1 in structure.sites:
        frac1 = numpy.array([site1.a, site1.b, site1.c])
        sites = structure.get_neighbors(site1, float(args[2]))
        nsite = len(sites)
        dcart = numpy.ndarray((nsite, 3), dtype=float)
        isite = 0
        for site2 in sites:
            frac2 = numpy.array([site2.a, site2.b, site2.c])
            dfrac = frac2[:] - frac1[:]
            dcart[isite, 0:3] = numpy.array(structure.lattice.get_cartesian_coords(dfrac))
            isite += 1

        ibase0 = 0
        sumdcart_max = 0.0
        for ibase in range(nbase):
            dcart_rot = numpy.dot(dcart, base[ibase, 0:3, 0:3])
            sumdcart = 0.0
            for isite in range(nsite):
                sumdcart += numpy.max(numpy.abs(dcart_rot[isite, 0:3]))
            if sumdcart > sumdcart_max:
                sumdcart_max = sumdcart
                ibase0 = ibase

        iaxis0 = 0
        dcart_rot = numpy.dot(dcart, base[ibase0, 0:3, 0:3])
        dcart_max = 0.0
        for isite in range(nsite):
            for iaxis in range(3):
                dcart_abs = abs(dcart_rot[isite, iaxis])
                if dcart_abs > dcart_max:
                    dcart_max = dcart_abs
                    iaxis0 = iaxis

        print(site1, nsite)
        print(site1.a, site1.b, site1.c, end=" ")
        for iaxis in range(3):
            iaxis1 = (iaxis0 + iaxis + 1) % 3
            dfrac = structure.lattice.get_fractional_coords(base[ibase0, 0:3, iaxis1])
            print(base[ibase0, 0, iaxis1], base[ibase0, 1, iaxis1], base[ibase0, 2, iaxis1], end=" ")
            print(icount*3+iaxis+1,
                  dfrac[0] * structure.lattice.a,
                  dfrac[1] * structure.lattice.b,
                  dfrac[2] * structure.lattice.c, 0, file=f)
            print(icount+1, 0, 0, 0, 0, file=f)
            print(0, 0, 0, 0, 0, file=f)
        print("")

        icount += 1
    print(0, 0, 0, 0, 0, file=f)

    print("VECTT", file=f)
    jcount = 1
    for site1 in structure.sites:
        sites = structure.get_neighbors(site1, 2.6)
        print(jcount,   0.5, 255,   0,   0, 0, file=f)
        print(jcount+1, 0.5,   0, 255,   0, 0, file=f)
        print(jcount+2, 0.5,   0,   0, 255, 0, file=f)
        jcount += 3
    print("0 0 0 0 00 0 0 0 0", file=f)
