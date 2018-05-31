#!/usr/bin/python3
import sys
import seekpath
import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy


if __name__ == '__main__':

    matcher = StructureMatcher()
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
        print(cif_file[0:len(cif_file) - 4], structure2.volume)
