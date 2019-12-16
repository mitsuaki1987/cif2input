#!/usr/bin/python3
import seekpath
import glob
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy
import re
import seekpath.hpkot


def main():
    with open("list.txt", "r") as f:
        files = f.readlines()
    matcher = StructureMatcher()
    #
    # Read All files specified as command-line arguments
    #
    for cif_file in files:
        cif_file = cif_file.strip("\n")
        print("Reading "+cif_file+" ... ", end="")
        #
        # PyMatGen structure from CIF file
        #
        try:
            structure = pymatgen.Structure.from_file(cif_file)
        except ValueError:
            print("Invalid structure.")
            continue
        except AssertionError:
            print("Invalid data.")
            continue
        #
        # Remove oxidation state. Excepting "D" (deuterium)
        #
        try:
            structure.remove_oxidation_states()
        except ValueError:
            print("Invalid element.")
            continue
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
        except seekpath.hpkot.SymmetryDetectionError:
            print("Except seekpath.hpkot.SymmetryDetection: Spglib could not detect the symmetry of the system")
            continue
        except ValueError:
            print("Problem creating primitive cell, I found the following group of atoms with len != 1: (0, 5)")
            continue
        #
        # Number of atoms
        #
        # nat_cut = 10
        # if structure2.num_sites > nat_cut:
        #     print("Number of atoms > " + str(nat_cut))
        #     continue
        #
        #
        #
        # dismat = structure2.distance_matrix
        # print(structure2.species)
        #
        # This structure is the same or not as the known structures
        #
        known = False
        for xsf_file in glob.glob("*.xsf"):
            known_structure = pymatgen.Structure.from_file(xsf_file)
            if matcher.fit(structure2, known_structure):
                print("Same as " + xsf_file)
                known = True
                break
        #
        # If it is new structure, save it as XSF
        #
        if not known:
            #
            xsf_file = ""
            for spc in structure2.composition.elements:
                xsf_file += str(spc) + "_" + str(structure2.species.count(spc))
            xsf_file += "-" + re.sub("\D", "", cif_file.split("/")[-1])+".xsf"
            print("Write to "+xsf_file)
            structure2.to(fmt="xsf", filename=xsf_file)


main()
