#!/usr/bin/python3
import sys
import glob
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy
import re


def main():
    matcher = StructureMatcher()
    args = sys.argv
    #
    # Read All files specified as command-line arguments
    #
    for ifile in range(len(args)-1):
        rx_file = args[ifile+1]
        print("Reading "+rx_file+" ... ", end="")
        #
        # PyMatGen structure from CIF file
        #
        with open(rx_file) as f:
            lines = f.readlines()
            if "Begin final coordinates\n" in lines:
                #
                # Parse from QE output
                #
                start = lines.index("Begin final coordinates\n")
                end = lines.index("End final coordinates\n")
                start_cell = lines[start:end].index("CELL_PARAMETERS (angstrom)\n") + start
                start_atom = lines[start:end].index("ATOMIC_POSITIONS (crystal)\n") + start
                #
                # Cell parameter
                #   avec[0,0:3] = a0
                #   avec[1,0:3] = a1
                #   avec[2,0:3] = a2
                #
                avec = numpy.array([lines[start_cell+1].strip().split(),
                                    lines[start_cell+2].strip().split(),
                                    lines[start_cell+3].strip().split()], numpy.float_)
                #
                natom = end - 1 - start_atom
                atom = numpy.zeros(natom, numpy.int_)
                pos = numpy.zeros((natom, 3), numpy.float_)
                iatom = 0
                for atom_pos in lines[start_atom+1:end]:
                    atom[iatom] = pymatgen.core.Element(atom_pos.strip().split()[0]).number
                    pos[iatom, 0:3] = numpy.array(atom_pos.strip().split()[1:4], numpy.float_)
                    iatom += 1
                #
                # PyMatGen structure
                #
                structure2 = pymatgen.core.Structure(avec, atom, pos)
                #
                # This structure is the same or not as the known structures
                #
                known = False
                for xsf_file in glob.glob("*.xsf"):
                    known_structure = pymatgen.core.Structure.from_file(xsf_file)
                    if matcher.fit(structure2, known_structure):
                        print("Same as " + xsf_file)
                        known = True
                        break
                #
                # If it is new structure, save it as XSF
                #
                if not known:
                    #
                    xsf_file = rx_file.replace("/", "")[:-6] + ".xsf"
                    print("Write to " + xsf_file)
                    structure2.to(fmt="xsf", filename=xsf_file)
            else:
                print("Not converged.")


main()
