#!/usr/bin/python3
import seekpath
import glob
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy
import re
import seekpath.hpkot
import CifFile
import sys


def main():
    #
    args = sys.argv
    with open(str(args[1]), "r") as f:
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
            structure = pymatgen.core.Structure.from_file(input_file)
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
        # Treatment for implicit hydrogen for CIF file
        #
        if input_file.split(".")[-1] == "cif":
            cf = CifFile.ReadCif(input_file)

            attached_hydrogens = cf.get_all("_atom_site_attached_hydrogens")
            if attached_hydrogens != ['0'] and attached_hydrogens != [] and attached_hydrogens != ['.']:
                print("Implicit hydrogen.", attached_hydrogens)
                continue

            try:
                formula_input = str(cf.get_all("_chemical_formula_sum")[0]).split(" ")
            except IndexError:
                print("No chemical formula.")
                continue
            formula_data = str(structure.composition.hill_formula).split(" ")
            try:
                formula_input_dict = {re.sub("[^a-zA-Z]", "", chem):
                                      1.0 if re.sub("[a-zA-Z]", "", chem) == "" else float(re.sub("[a-zA-Z]", "", chem))
                                      for chem in formula_input}
            except ValueError:
                print("Invalid chemical formula string.")
                continue
            formula_data_dict = {re.sub("[^a-zA-Z]", "", chem):
                                 1.0 if re.sub("[a-zA-Z]", "", chem) == "" else float(re.sub("[a-zA-Z]", "", chem))
                                 for chem in formula_data}
            formula_input = re.sub(" ", "", str(cf.get_all("_chemical_formula_sum")[0]))
            formula_data = re.sub(" ", "", str(structure.composition.hill_formula))
            if len(formula_data_dict) != len(formula_input_dict):
                print("Number of elements are invalid. ", formula_input, formula_data)
                continue
            else:
                elem_switch = -1
                ratio0 = 0.0
                for elem in formula_input_dict:
                    if elem in formula_data_dict:
                        ratio = formula_input_dict[elem] / formula_data_dict[elem]
                        if elem_switch == -1:
                            ratio0 = ratio
                            elem_switch = 0
                        else:
                            if abs(ratio - ratio0) > 0.001:
                                elem_switch = 1
                                break
                    else:
                        elem_switch = 2
                        break
                if elem_switch == 1:
                    print("Chemical formula does not match.", formula_input, formula_data)
                    continue
                elif elem_switch == 2:
                    print("Elements does not match", formula_input, formula_data)
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
                                               [pymatgen.core.Element(str(spc)).number for spc in structure.species]))
            structure2 = pymatgen.core.Structure(skp["primitive_lattice"],
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
        output_file = re.sub(" ", "", str(structure2.composition.alphabetical_formula))
        #
        # This structure is the same or not as the known structures
        #
        known = False
        print(glob.glob(output_file + "-*.xsf"), " ... ", end="")
        for known_file in glob.glob(output_file + "-*.xsf"):
            known_structure = pymatgen.core.Structure.from_file(known_file)
            if matcher.fit(structure2, known_structure):
                print("Same as " + known_file)
                known = True
                break
        #
        # If it is new structure, save it as XSF
        #
        if not known:
            #
            if input_file.split(".")[-1] == "xsf":
                output_file += "-" + re.sub("[^0-9]", "", input_file.split("-")[-1]) + ".xsf"
            else:
                output_file += "-" + re.sub("[^0-9]", "", input_file.split("/")[-1])+".xsf"
            print("Write to "+output_file)
            structure2.to(fmt="xsf", filename=output_file)


main()
