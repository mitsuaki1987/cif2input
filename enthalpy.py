#!/usr/bin/python3
import sys

if __name__ == '__main__':

    args = sys.argv

    if args[1] == "sg15":
        from sg15 import psen_dict
    elif args[1] == "pslibrary":
        from pslibrary import psen_dict
    else:
        from sssp import psen_dict

    for ifile in range(len(args)-2):
        rx_file = args[ifile+2]

        with open(rx_file) as f:
            lines = f.readlines()
            if "Begin final coordinates\n" in lines:
                start = lines.index("Begin final coordinates\n")
                end = lines.index("End final coordinates\n")
                natom = len(lines[start+10:end])
                enthalpy = float(lines[start-1].strip().split()[3]) / float(natom)
                #
                # isolated
                #
                enthalpy0 = 0.0
                iatom = 0
                for atom_pos in lines[start+10:end]:
                    atom = str(atom_pos.strip().split()[0])
                    enthalpy0 += psen_dict[atom]
                    iatom += 1
                enthalpy0 /= float(natom)

                print(rx_file, enthalpy, enthalpy0, enthalpy - enthalpy0)
