#!/usr/bin/python3
import sys
from sssp import atomwfc_dict


def main():
    args = sys.argv
    with open(str(args[1]), "r") as f:
        name = f.readlines()
    print("{")

    nname = len(name)

    for iname in name:

        iname = iname.strip("\n")
        print("  \"" + iname + "\" : {")

        with open("scfout/scf_" + iname + ".out", "r") as f:
            output_lines = f.readlines()
        etot = 0.0
        magt = 0.0
        nat = 1
        for output_line in output_lines:
            nat_words = output_line.split("number of atoms/cell      =")
            etot_words = output_line.split("!    total energy              =")
            magt_words = output_line.split("total magnetization       =")
            if len(nat_words) > 1:
                nat = int(nat_words[1].split()[0])
            elif len(etot_words) > 1:
                etot = float(etot_words[1].split()[0])
            elif len(magt_words) > 1:
                magt = float(magt_words[1].split()[0])

        with open("isoen/" + iname + "_iso.dat", "r") as f:
            eiso = float(f.readline())

        print("    \"eform\" : " + str((etot-eiso)/float(nat)) + ",")
        print("    \"magt\" : " + str(magt/nat) + ",")

        with open("dos/" + iname + ".pdos_tot", "r") as f:
            output_lines = f.readlines()[1].split()
        dosf = float(output_lines[1]) + float(output_lines[2])
        print("    \"dosf\" : " + str(dosf/float(nat)) + ",")

        with open("scfin/scf_" + iname + ".in", "r") as f:
            output_lines = f.readlines()
        types = []
        reading = False
        for output_line in output_lines:
            if output_line == "ATOMIC_SPECIES\n":
                reading = True
            elif output_line == "ATOMIC_POSITIONS crystal\n":
                break
            elif reading:
                types.append(output_line.split()[0])
        atoms = []
        reading = False
        for output_line in output_lines:
            if output_line == "ATOMIC_POSITIONS crystal\n":
                reading = True
            elif output_line == "K_POINTS automatic\n":
                break
            elif reading:
                atoms.append(output_line.split()[0])

        ntype = len(types)

        pdos = {}
        for iat in range(len(atoms)):
            if atoms[iat] not in pdos:
                pdos[atoms[iat]] = {}
            for iwfc in range(len(atomwfc_dict[atoms[iat]])):
                pnum = atomwfc_dict[atoms[iat]][iwfc][0][0]
                ang = atomwfc_dict[atoms[iat]][iwfc][0][1]
                if pnum not in pdos[atoms[iat]]:
                    pdos[atoms[iat]][pnum] = {}
                if ang not in pdos[atoms[iat]][pnum]:
                    pdos[atoms[iat]][pnum][ang] = []
                with open("pdos_atm/" + iname + ".pdos_atm#" + str(iat + 1) + "(" + atoms[iat] + ")_wfc#"
                          + str(iwfc + 1) + "(" + ang.lower() + ")", "r") as f:
                    output_lines = f.readlines()[1].split()
                pdos[atoms[iat]][pnum][ang].append(float(output_lines[1]) + float(output_lines[2]))

        print("    \"type\" : {")
        for itype in types:
            print("      \"" + itype + "\" : {")
            print("        \"nat\" : " + str(atoms.count(itype)) + ",")

            npnam = len(pdos[itype])
            for ipnam in pdos[itype].keys():
                print("        \"" + ipnam + "\" : {")
                nang = len(pdos[itype][ipnam])
                for iang in pdos[itype][ipnam].keys():
                    print("          \"" + iang + "\" : {")
                    print("            \"pdos\" : " + str(pdos[itype][ipnam][iang]))

                    nang += -1
                    if nang == 0:
                        print("          }")
                    else:
                        print("          },")

                npnam += -1
                if npnam == 0:
                    print("        }")
                else:
                    print("        },")
            ntype += -1
            if ntype == 0:
                print("      }")
            else:
                print("      },")
        print("    }")

        nname += -1
        if nname == 0:
            print("  }")
        else:
            print("  },")

    print("}")


main()
