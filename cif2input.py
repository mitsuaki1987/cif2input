#!/usr/bin/python3
import sys
import pymatgen
from structure2input import structure2input


def main():
    args = sys.argv
    if len(args) < 2:
        print("Usage:")
        print("$ cif2input.py cif-file [pseudo_kind] [host] [rel] [uc] [dk_path] [dq_grid]")
        print("Default:")
        print("$ cif2input.py cif-file sg15 enaga 0 1 0.1 0.3359385398275")
        print("pseudo_kind : sg15, pslibrary, sssp, ssspsol")
        print("host : enaga, ohtaka, wisteria")
        print("rel : 0 (scaler/non relativistic) or 1 (spin-orbit)")
        print("uc : 0 (do not) or 1 (reduce to primitive unit-cell)")
        exit(0)
    #
    # CIF parser
    #
    structure = pymatgen.core.Structure.from_file(args[1])
    #
    # Default value
    #
    dk_path = 0.1
    dq_grid = 0.27
    pseudo_kind = "sssp"
    host = "wisteria"
    irel = 0
    iuc = 1
    #
    if len(args) > 2:
        pseudo_kind = args[2]
        if len(args) > 3:
            host = args[3]
            if len(args) > 4:
                irel = int(args[4])
                if len(args) > 5:
                    iuc = int(args[5])
                    if len(args) > 6:
                        dk_path = float(args[6])
                        if len(args) > 7:
                            dq_grid = float(args[7])
    rel = False
    if irel == 0:
        rel = False
    elif irel == 1:
        rel = True
    else:
        print("Invalid rel :", irel)
        exit(irel)
    #
    uc = False
    if iuc == 0:
        uc = False
    elif iuc == 1:
        uc = True
    else:
        print("Invalid uc :", irel)
        exit(irel)
    #
    print("  dk for band : {0}".format(dk_path))
    print("  dq for grid : {0}".format(dq_grid))
    print("  Pseudo kind : ", pseudo_kind)
    print("         Host : ", host)
    print(" Relativistic : ", rel)
    print(" Reduce to UC : ", uc)

    structure.remove_oxidation_states()

    structure2input(structure, dk_path, dq_grid, pseudo_kind, host, rel, uc)


main()
