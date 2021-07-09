#!/usr/bin/python3
import sys
import seekpath
import pymatgen
import numpy
from pymatgen.core.periodic_table import get_el_sp
import math


def main():
    e_ion_dict ={
        "H": 1311.3,
        "He": 2361.3,
        "Li": 519.9,
        "Be": 898.8,
        "B": 800.2,
        "C": 1085.7,
        "N": 1401.5,
        "O": 1313.1,
        "F": 1680.0,
        "Ne": 2079.4,
        "Na": 495.6,
        "Mg": 737.3,
        "Al": 577.5,
        "Si": 786.0,
        "P": 1011.2,
        "S": 999.0,
        "Cl": 1254.9,
        "Ar": 1519.6,
        "K": 418.5,
        "Ca": 589.4,
        "Sc": 630.8,
        "Ti": 657.8,
        "V": 650.1,
        "Cr": 652.4,
        "Mn": 716.8,
        "Fe": 759.1,
        "Co": 758.1,
        "Ni": 736.2,
        "Cu": 745.0,
        "Zn": 905.8,
        "Ga": 578.7,
        "Ge": 760.0,
        "As": 946.2,
        "Se": 940.4,
        "Br": 1142.0,
        "Kr": 1350.0,
        "Rb": 402.8,
        "Sr": 549.0,
        "Y": 615.4,
        "Zr": 659.7,
        "Nb": 663.6,
        "Mo": 684.8,
        "Tc": 702.2,
        "Ru": 710.3,
        "Rh": 719.5,
        "Pd": 803.5,
        "Ag": 730.5,
        "Cd": 867.5,
        "In": 558.0,
        "Sn": 708.2,
        "Sb": 833.3,
        "Te": 869.0,
        "I": 1008.3,
        "Xe": 1170.0,
        "Cs": 375.5,
        "Ba": 502.5,
        "La": 541.1,
        "Ce": 540.1,
        "Pr": 526.6,
        "Nd": 531.5,
        "Pm": 536,
        "Sm": 540.1,
        "Eu": 546.9,
        "Gd": 594.2,
        "Tb": 569,
        "Dy": 567,
        "Ho": 574,
        "Er": 581,
        "Tm": 589,
        "Yb": 603,
        "Lu": 513,
        "Hf": 575.2,
        "Ta": 760.1,
        "W": 769.7,
        "Re": 759.1,
        "Os": 819.8,
        "Ir": 868.1,
        "Pt": 868.1,
        "Au": 889.3,
        "Hg": 1006.0,
        "Tl": 588.9,
        "Pb": 715.2,
        "Bi": 702.9,
        "Po": 813.1,
        "At": 916.3,
        "Rn": 1036.5,
        "Fr": 380,
        "Ra": 509.3,
        "Ac": 665.5,
        "Th": 670.4,
        "Pa": 0.0,
        "U": 686.4,
        "Np": 0.0,
        "Pu": 538.3,
        "Am": 0.0,
        "Cm": 581,
        "Bk": 600,
        "Cf": 610,
        "Es": 619,
        "Fm": 628.5
    }

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
        structure = pymatgen.core.Structure.from_file(input_file)
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
                                               [pymatgen.core.Element(str(spc)).number for spc in structure.species]))
            structure2 = pymatgen.core.Structure(skp["primitive_lattice"],
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

        desc = []
        desc.append(vol)  # 2
        desc.append(ave_period)  # 3
        desc.append(ave_num_d)  # 4
        desc.append(ave_oxidization)  # 5
        desc.append(var_period)  # 6
        desc.append(var_num_d)  # 7
        desc.append(var_oxidization)  # 8
        desc.append(distance_min)  # 9
        
        print(input_file, end=" ")
        for i_desc in range(len(desc)):
            print(desc[i_desc], end=" ")
        for i_desc in range(len(desc)):
            for j_desc in range(i_desc+1):
                print(desc[i_desc]*desc[j_desc], end=" ")
        print("")

main()
