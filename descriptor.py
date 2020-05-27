#!/usr/bin/python3
import sys
import seekpath
import pymatgen
import numpy
from pymatgen.core.periodic_table import get_el_sp
import math


def main():
    valence_dict = {
        "H":  [1,  0,  0, 0],
        "He": [1,  0,  0, 0],
        "Li": [2,  1,  0, 0],
        "Be": [2,  2,  0, 0],
        "B":  [2,  3,  0, 0],
        "C":  [2,  4,  0, 0],
        "N":  [2, -3,  0, 0],
        "O":  [2, -2,  0, 0],
        "F":  [2, -1,  0, 0],
        "Ne": [2,  0,  0, 0],
        "Na": [3,  1,  0, 0],
        "Mg": [3,  2,  0, 0],
        "Al": [3,  3,  0, 0],
        "Si": [3,  4,  0, 0],
        "P":  [3, -3,  0, 0],
        "S":  [3, -2,  0, 0],
        "Cl": [3, -1,  0, 0],
        "Ar": [3,  0,  0, 0],
        "K":  [4,  1,  0, 0],
        "Ca": [4,  2,  0, 0],
        "Sc": [4,  2,  1, 0],
        "Ti": [4,  2,  2, 0],
        "V":  [4,  2,  3, 0],
        "Cr": [4,  2,  4, 0],
        "Mn": [4,  2,  5, 0],
        "Fe": [4,  2,  6, 0],
        "Co": [4,  2,  7, 0],
        "Ni": [4,  2,  8, 0],
        "Cu": [4,  2,  9, 0],
        "Zn": [4,  2, 10, 0],
        "Ga": [4,  3,  0, 0],
        "Ge": [4,  4,  0, 0],
        "As": [4, -3,  0, 0],
        "Se": [4, -2,  0, 0],
        "Br": [4, -1,  0, 0],
        "Kr": [4,  0,  0, 0],
        "Rb": [5,  1,  0, 0],
        "Sr": [5,  2,  0, 0],
        "Y":  [5,  2,  1, 0],
        "Zr": [5,  2,  2, 0],
        "Nb": [5,  2,  3, 0],
        "Mo": [5,  2,  4, 0],
        "Tc": [5,  2,  5, 0],
        "Ru": [5,  2,  6, 0],
        "Rh": [5,  2,  7, 0],
        "Pd": [5,  2,  8, 0],
        "Ag": [5,  2,  9, 0],
        "Cd": [5,  2, 10, 0],
        "In": [5,  3,  0, 0],
        "Sn": [5,  4,  0, 0],
        "Sb": [5, -3,  0, 0],
        "Te": [5, -2,  0, 0],
        "I":  [5, -1,  0, 0],
        "Xe": [5,  0,  0, 0],
        "Cs": [6,  1,  0, 0],
        "Ba": [6,  2,  0, 0],
        "La": [6,  2,  1, 0],
        "Ce": [6,  2,  1, 1],
        "Pr": [6,  2,  0, 3],
        "Nd": [6,  2,  0, 4],
        "Pm": [6,  2,  0, 5],
        "Sm": [6,  2,  0, 6],
        "Eu": [6,  2,  0, 7],
        "Gd": [6,  2,  0, 8],
        "Tb": [6,  2,  0, 9],
        "Dy": [6,  2,  0, 10],
        "Ho": [6,  2,  0, 11],
        "Er": [6,  2,  0, 12],
        "Tm": [6,  2,  0, 13],
        "Yb": [6,  2,  0, 14],
        "Lu": [6,  2,  1, 0],
        "Hf": [6,  2,  2, 0],
        "Ta": [6,  2,  3, 0],
        "W":  [6,  2,  4, 0],
        "Re": [6,  2,  5, 0],
        "Os": [6,  2,  6, 0],
        "Ir": [6,  2,  7, 0],
        "Pt": [6,  2,  8, 0],
        "Au": [6,  2,  9, 0],
        "Hg": [6,  2, 10, 0],
        "Tl": [6,  3,  0, 0],
        "Pb": [6,  4,  0, 0],
        "Bi": [6, -3,  0, 0],
        "Po": [6, -2,  0, 0],
        "At": [6, -1,  0, 0],
        "Rn": [6,  0,  0, 0],
        "Fr": [7,  1,  0, 0],
        "Ra": [7,  2,  0, 0],
        "Ac": [7,  2,  1, 0],
        "Th": [7,  2,  1, 1],
        "Pa": [7,  2,  0, 3],
        "U":  [7,  2,  0, 4],
        "Np": [7,  2,  0, 5],
        "Pu": [7,  2,  0, 6],
        "Am": [7,  2,  0, 7],
        "Cm": [7,  2,  0, 8],
        "Bk": [7,  2,  0, 9],
        "Cf": [7,  2,  0, 10],
        "Es": [7,  2,  0, 11],
        "Fm": [7,  2,  0, 12],
    }

    args = sys.argv
    with open(str(args[1]), "r") as f:
        input_list = f.readlines()
    #
    # Read All files specified as command-line arguments
    #
    for input_file in input_list:
        input_file = input_file.strip("\n")
        #
        # PyMatGen structure from CIF file
        #
        structure = pymatgen.Structure.from_file(input_file)
        structure.remove_oxidation_states()
        #
        # Refine 3-folded Wyckoff position
        #
        frac_coord2 = numpy.array(structure.frac_coords)
        for ipos in range(len(frac_coord2)):
            for iaxis in range(3):
                coord3 = frac_coord2[ipos, iaxis] * 6.0
                if abs(round(coord3) - coord3) < 0.001:
                    frac_coord2[ipos, iaxis] = float(round(coord3)) / 6.0
        #
        # Disordered structure raises AttributeError. So it is skipped.
        #
        try:
            skp = seekpath.get_explicit_k_path((structure.lattice.matrix, frac_coord2,
                                               [pymatgen.Element(str(spc)).number for spc in structure.species]))
            structure2 = pymatgen.Structure(skp["primitive_lattice"],
                                            skp["primitive_types"], skp["primitive_positions"])
        except AttributeError:
            print("Fractional occupancy, may be disordered.")
            continue
        #
        nat = len(structure2.atomic_numbers)
        atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
        #
        # vol : Volume per atom
        #
        vol = structure2.volume / float(nat)
        #
        # Average : oxidization, d-electron, period
        #
        ave_period = 0.0
        ave_num_d = 0.0
        ave_oxidization = 0.0
        for iat in atom:
            ave_period += valence_dict[iat][0]
            ave_num_d += valence_dict[iat][2]
            ave_oxidization += valence_dict[iat][1] \
                + valence_dict[iat][2] \
                + valence_dict[iat][3]
        ave_period /= float(nat)
        ave_num_d /= float(nat)
        ave_oxidization /= float(nat)
        #
        # Variance : oxidization, d-electron, period
        #
        var_period = 0.0
        var_num_d = 0.0
        var_oxidization = 0.0
        for iat in atom:
            var_period += (ave_period - valence_dict[iat][0])**2
            var_num_d += (ave_num_d - valence_dict[iat][2])**2
            var_oxidization += (ave_oxidization
                                - valence_dict[iat][1]
                                - valence_dict[iat][2]
                                - valence_dict[iat][3])**2
        var_period = math.sqrt(var_period / float(nat))
        var_num_d = math.sqrt(var_num_d / float(nat))
        var_oxidization = math.sqrt(var_oxidization / float(nat))
        #
        # distance_min : Minimum distance
        #
        if nat == 1:
            distance_min = numpy.min(structure2.lattice.abc)
        else:
            distance_matrix = structure2.distance_matrix
            distance_min = distance_matrix[0, 1]
            for iat in range(nat):
                for jat in range(iat+1, nat):
                    distance_min = min(distance_min, distance_matrix[iat, jat])
        
        print(input_file,
              vol, ave_period, ave_num_d, ave_oxidization, var_period, var_num_d, var_oxidization, distance_min,
              vol**2, ave_period**2, ave_num_d**2, ave_oxidization**2,
              var_period**2, var_num_d**2, var_oxidization**2, distance_min**2,
              1.0/vol, 1.0/ave_period, 1.0/distance_min,
              1.0/vol**2, 1.0/ave_period**2, 1.0/distance_min**2
              )


main()
