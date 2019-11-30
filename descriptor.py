#!/usr/bin/python3
import sys
import seekpath
import pymatgen
import numpy
from pymatgen.core.periodic_table import get_el_sp


def main():
    valence_dict = {
        "H": 0,
        "He": 0,
        "Li": 1,
        "Be": 1.5,
        "B": 3,
        "C": 0,
        "N": 1,
        "O": -0.125,
        "F": -0.5,
        "Ne": 0,
        "Na": 1,
        "Mg": 2,
        "Al": 3,
        "Si": 0.5,
        "P": 0.833333,
        "S": 1.42857,
        "Cl": 3.125,
        "Ar": 0,
        "K": 1,
        "Ca": 2,
        "Sc": 3,
        "Ti": 3,
        "V": 2.8,
        "Cr": 2.75,
        "Mn": 3,
        "Fe": 2.75,
        "Co": 1,
        "Ni": 1.66667,
        "Cu": 1.5,
        "Zn": 2,
        "Ga": 3,
        "Ge": 3,
        "As": 1.66667,
        "Se": 2.66667,
        "Br": 2.5,
        "Kr": 2,
        "Rb": 1,
        "Sr": 2,
        "Y": 3,
        "Zr": 2,
        "Nb": 3,
        "Mo": 3.33333,
        "Tc": 3,
        "Ru": 3,
        "Rh": 2.5,
        "Pd": 3,
        "Ag": 1.5,
        "Cd": 2,
        "In": 2,
        "Sn": 3,
        "Sb": 1.66667,
        "Te": 2.5,
        "I": 2.5,
        "Xe": 4,
        "Cs": 1,
        "Ba": 2,
        "La": 3,
        "Ce": 3.5,
        "Pr": 3.5,
        "Nd": 3,
        "Pm": 3,
        "Sm": 2.5,
        "Eu": 2.5,
        "Gd": 3,
        "Tb": 3.5,
        "Dy": 3,
        "Ho": 3,
        "Er": 3,
        "Tm": 2.5,
        "Yb": 2.5,
        "Lu": 3,
        "Hf": 4,
        "Ta": 5,
        "W": 3.33333,
        "Re": 3.71429,
        "Os": 3,
        "Ir": 2.14286,
        "Pt": 2,
        "Au": 2,
        "Hg": 1.5,
        "Tl": 2,
        "Pb": 2,
        "Bi": 4,
        "Po": 2.5,
        "At": 3,
        "Rn": 5,
        "Fr": 1.5,
        "Ra": 2,
        "Ac": 3,
        "Th": 4,
        "Pa": 4.5,
        "U": 4.5,
        "Np": 5,
        "Pu": 4.5,
        "Am": 4.5,
        "Cm": 3.5,
        "Bk": 3.5,
        "Cf": 3,
        "Es": 3,
        "Fm": 3,
    }

    args = sys.argv
    #
    # Read All files specified as command-line arguments
    #
    for ifile in range(len(args)-1):
        cif_file = args[ifile+1]
        #
        # PyMatGen structure from CIF file
        #
        structure = pymatgen.Structure.from_file(cif_file)
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
        # density : 1 / Volume
        #
        density = 1.0 / vol
        #
        # oxidization : Average oxidization
        #
        oxidization = 0
        for iat in atom:
            oxidization += valence_dict[iat]
        oxidization /= float(nat)
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

        print(cif_file, vol, density, oxidization, distance_min)


main()
