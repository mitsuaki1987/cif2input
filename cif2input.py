#!/usr/bin/python3
import sys
import pymatgen
from structure2input import structure2input

if __name__ == '__main__':

    args = sys.argv
    if len(args) < 2:
        print("Usage:")
        print("$ cif2input.py cif-file [dk_path] [dq_grid] [pseudo_kind] [queue] [rel]")
        print("Default:")
        print("$ cif2input.py cif-file 0.1 0.3359385398275 sg15 F4cpue")
        exit(0)
    #
    # CIF parser
    #
    structure = pymatgen.Structure.from_file(args[1])
    #
    # Default value
    #
    dk_path = 0.1
    dq_grid = 0.3359385398275
    pseudo_kind = "sg15"
    pseudo_dir = "/work/i0012/i001200/pseudo/"
    queue = "F4cpue"
    rel = False
    #
    if len(args) > 2:
        dk_path = float(args[2])
        if len(args) > 3:
            dq_grid = float(args[3])
            if len(args) > 4:
                pseudo_kind = args[4]
                if len(args) > 5:
                    queue = args[5]
                    if len(args) > 6:
                        rel = True
    #
    print("  dk for band : {0}".format(dk_path))
    print("  dq for grid : {0}".format(dq_grid))
    print("  Pseudo kind is ", pseudo_kind)
    print("  Pseudo is at ", pseudo_dir)

    structure.remove_oxidation_states()

    structure2input(structure, dk_path, dq_grid, pseudo_kind, pseudo_dir, queue, rel)
