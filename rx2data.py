#!/usr/bin/python3
import sys
import glob
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy


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
                avec = numpy.array([lines[start+5].strip().split(),
                                    lines[start+6].strip().split(),
                                    lines[start+7].strip().split()], numpy.float_)
                natom = len(lines[start+10:end])
                atom = numpy.zeros(natom, numpy.int_)
                pos = numpy.zeros((natom, 3), numpy.float_)
                iatom = 0
                for atom_pos in lines[start+10:end]:
                    atom[iatom] = pymatgen.Element(atom_pos.strip().split()[0]).number
                    pos[iatom, 0:3] = numpy.array(atom_pos.strip().split()[1:4], numpy.float_)
                    iatom += 1
                #
                # PyMatGen structure
                #
                structure2 = pymatgen.Structure(avec, atom, pos)
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
                    xsf_file = rx_file[0:len(rx_file) - 9] + ".xsf"
                    print("Write to " + xsf_file)
                    structure2.to(fmt="xsf", filename=xsf_file)
            else:
                print("Not converged.")
