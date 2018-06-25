#!/usr/bin/python3
import pslibrary
import sg15
import sssp

ii = 0
for lib in pslibrary, sg15, sssp:
    ii += 1
    for atom in lib.pseudo_dict.keys():
        with open(str(ii)+"/"+atom+".in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'scf\'", file=f)
            print("  pseudo_dir = \'/work/i0012/i001200/pseudo/\'", file=f)
            print("      prefix = \'%s\'" % atom, file=f)
            print("/", file=f)
            print("&SYSTEM", file=f)
            print("       ibrav = 1", file=f)
            print("   celldm(1) = 20.0", file=f)
            print("         nat = 1", file=f)
            print("        ntyp = 1", file=f)
            print("     ecutwfc = %f" % lib.ecutwfc_dict[atom], file=f)
            print("     ecutrho = %f" % lib.ecutrho_dict[atom], file=f)
            print(" occupations = \'smearing\'", file=f)
            print("    smearing = \'m-p\'", file=f)
            print("     degauss = 0.05", file=f)
            print("    assume_isolated = \'m-t\'", file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" conv_thr = 1.0e-8", file=f)
            print(" mixing_beta = 0.3", file=f)
            print("/", file=f)
            print("ATOMIC_SPECIES", file=f)
            print(" %s 1.0 %s" % (atom, lib.pseudo_dict[atom]), file=f)
            print("ATOMIC_POSITIONS crystal", file=f)
            print(" %s 0.0 0.0 0.0" % atom, file=f)
            print("K_POINTS gamma", file=f)
