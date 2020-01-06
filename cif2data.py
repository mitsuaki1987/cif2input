#!/usr/bin/python3
import seekpath
import glob
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy
import re
import seekpath.hpkot


def main():
    #
    # Collect known structures
    #
    known_structures = []
    for known_file in glob.glob("*.xsf"):
        known_structures.append(pymatgen.Structure.from_file(known_file))
    print("Number of known structure : ", len(known_structures))
    #
    with open("list.txt", "r") as f:
        files = f.readlines()
    matcher = StructureMatcher()
    #
    # Read All files specified as command-line arguments
    #
    for input_file in files:
        input_file = input_file.strip("\n")
        print("Reading "+input_file+" ... ", end="")
        #
        # PyMatGen structure from CIF file
        #
        try:
            structure = pymatgen.Structure.from_file(input_file)
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
            if len(skp["primitive_types"]) != 1:
                structure2.merge_sites(tol=0.01, mode="average")
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
        # This structure is the same or not as the known structures
        #
        known = False
        for known_structure in known_structures:
            if matcher.fit(structure2, known_structure):
                print("Same as known structure")
                known = True
                break
        #
        # If it is new structure, save it as XSF
        #
        if not known:
            #
            output_file = re.sub(" ", "", str(structure2.composition.alphabetical_formula))
            #
            if input_file.split(".")[-1] == "xsf":
                output_file += "-" + re.sub("\D", "", input_file.split("-")[-1]) + ".xsf"
            else:
                output_file += "-" + re.sub("\D", "", input_file.split("/")[-1])+".xsf"
            print("Write to "+output_file)
            structure2.to(fmt="xsf", filename=output_file)
            known_structures.append(structure2)


main()
