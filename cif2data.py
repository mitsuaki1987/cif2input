#!/usr/bin/python3
import sys
import seekpath
import glob
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher


if __name__ == '__main__':

    matcher = StructureMatcher()
    args = sys.argv
    for ifile in range(len(args)-1):
        cif_file = args[ifile+1]
        print("Examining "+cif_file+"... ", end="")
        structure = pymatgen.Structure.from_file(cif_file)
        structure.remove_oxidation_states()
        skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                                            [pymatgen.Element(str(spc)).number for spc in structure.species]))
        structure2 = pymatgen.Structure(skp["primitive_lattice"], skp["primitive_types"], skp["primitive_positions"])

        known = False
        for xsf_file in glob.glob("*.xsf"):
            known_structure = pymatgen.Structure.from_file(xsf_file)
            if matcher.fit(structure2, known_structure):
                print(" Same as " + xsf_file)
                known = True
                break

        if not known:
            xsf_file = cif_file[0:len(cif_file) - 4]+".xsf"
            print(" Write to "+xsf_file)
            structure2.to(fmt="xsf", filename=xsf_file)
