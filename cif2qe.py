#!/usr/bin/python3
import sys
import pymatgen
from sssp import pseudo_dict

args = sys.argv

structure = pymatgen.Structure.from_file(args[1])

structure = structure.get_primitive_structure()
structure = structure.get_primitive_structure()
structure = structure.get_primitive_structure()
structure = structure.get_primitive_structure()

structure.remove_oxidation_states()

avec = structure.lattice.matrix
nat = len(structure.species)
typ = set(structure.species)
ntyp = len(typ)

print(structure)

with open("scf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'scf\'", file=f)
    print("      outdir = \'./\'", file=f)
    print("  pseudo_dir = \'../pseudo/\'", file=f)
    print("      prefix = \'%s\'" % args[2], file=f)
    print("/", file=f)

    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %d" % 30.0, file=f)
    print("     ecutrho = %d" % 300.0, file=f)
    print(" occupations = \'tetrahedra_opt\'", file=f)
    print("/", file=f)

    print("&ELECTRONS", file=f)
    print("/", file=f)

    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii][0],avec[ii][1],avec[ii][2]), file=f)

    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)

    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            structure.species[iat],
            structure.frac_coords[iat][0],structure.frac_coords[iat][1],structure.frac_coords[iat][2]), file=f)

    print("K_POINTS automatic", file=f)
    print(" 4 4 4 0 0 0", file=f)
