#!/usr/bin/python3
import sys
import pymatgen
from structure2input import structure2input

if __name__ == '__main__':

    args = sys.argv
    if len(args) < 2:
        print("Usage:")
        print("$ cif2input cif-file [prefix] [dk_path] [dq_grid] [pseudo_kind] [pseudo_path]")
        print("Default:")
        print("$ cif2input cif-file cif-file 0.1 0.3359385398275 sssp ../pseudo/")
        exit(0)
    #
    # CIF parser
    #
    structure = pymatgen.Structure.from_file(args[1])
    #
    # Default value
    #
    prefix = args[1][0:len(args[1]) - 4]
    dk_path = 0.1
    dq_grid = 0.3359385398275
    pseudo_kind = "sssp"
    pseudo_dir = "/work/i0012/i001200/pseudo/"
    #
    if len(args) > 2:
        prefix = args[2]
        if len(args) > 3:
            dk_path = float(args[3])
            if len(args) > 4:
                dq_grid = float(args[4])
                if len(args) > 5:
                    pseudo_kind = args[5]
                    if len(args) > 6:
                        pseudo_dir = args[6]
    #
    print("  prefix : {0}".format(prefix))
    print("  dk for band : {0}".format(dk_path))
    print("  dq for grid : {0}".format(dq_grid))
    print("  Pseudo kind is ", pseudo_kind)
    print("  Pseudo is at ", pseudo_dir)

    structure.remove_oxidation_states()

    structure2input(structure, prefix, dk_path, dq_grid, pseudo_kind, pseudo_dir)
