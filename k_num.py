#!/usr/bin/python3
import pymatgen
from pymatgen.core.periodic_table import get_el_sp
from sssp import ecutwfc_dict, ecutrho_dict, valence_dict
import numpy
import sys
import math
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def main():
    #
    args = sys.argv
    with open(str(args[1]), "r") as f:
        input_list = f.readlines()
    #
    # Read previous result
    #
    for input_file in input_list:
        #
        input_file = input_file.strip("\n")
        #
        structure = pymatgen.core.Structure.from_file(input_file)
        avec = structure.lattice.matrix
        bvec = structure.lattice.reciprocal_lattice.matrix
        atom = [str(get_el_sp(iat)) for iat in structure.atomic_numbers]
        typ = sorted(set(atom))
        bz_volume = abs(- bvec[0][2]*bvec[1][1]*bvec[2][0]
                        + bvec[0][1]*bvec[1][2]*bvec[2][0]
                        + bvec[0][2]*bvec[1][0]*bvec[2][1]
                        - bvec[0][0]*bvec[1][2]*bvec[2][1]
                        - bvec[0][1]*bvec[1][0]*bvec[2][2]
                        + bvec[0][0]*bvec[1][1]*bvec[2][2])
        #
        # WFC and Rho cutoff
        #
        ecutwfc = 0.0
        ecutrho = 0.0
        for ityp in typ:
            if ecutwfc < ecutwfc_dict[str(ityp)]:
                ecutwfc = ecutwfc_dict[str(ityp)]
            if ecutrho < ecutrho_dict[str(ityp)]:
                ecutrho = ecutrho_dict[str(ityp)]
        numpw = 4.0 * math.pi * math.sqrt(ecutwfc) ** 3 / 3.0 / bz_volume
        #
        # k and q grid
        #  the number of grid proportional to the Height of b
        #  b_i * a_i / |a_i| = 2pi / |a_i| (a_i is perpendicular to other b's)
        #
        nk = numpy.zeros(3, numpy.int_)
        for ii in range(3):
            norm = numpy.sqrt(numpy.dot(avec[ii][:], avec[ii][:]))
            nk[ii] = round(2.0 * numpy.pi / norm / 0.15)
            if nk[ii] == 0:
                nk[ii] = 1
        spg_analysis = SpacegroupAnalyzer(structure)
        coarse = spg_analysis.get_ir_reciprocal_mesh(mesh=(nk[0], nk[1], nk[2]), is_shift=(0, 0, 0))
        dense = spg_analysis.get_ir_reciprocal_mesh(mesh=(nk[0]*2, nk[1]*2, nk[2]*2), is_shift=(0, 0, 0))
        #
        # Number of electrons
        #
        nbnd = 0
        for iat in atom:
            nbnd += valence_dict[iat]
        nbnd = nbnd / 2
        #
        # Summary
        #
        print(input_file, len(coarse), len(dense), numpw, nbnd)


main()
