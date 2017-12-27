#!/usr/bin/python3
import sys
import numpy
import pymatgen
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict

args = sys.argv

structure = pymatgen.Structure.from_file(args[1])

structure = structure.get_primitive_structure()
structure = structure.get_primitive_structure()
structure = structure.get_primitive_structure()
structure = structure.get_primitive_structure()

structure.remove_oxidation_states()

avec = structure.lattice.matrix
bvec = structure.lattice.reciprocal_lattice.matrix
print(bvec)
nat = len(structure.species)
typ = set(structure.species)
ntyp = len(typ)
ecutwfc = 0.0
ecutrho = 0.0
for ityp in typ:
    if ecutwfc < ecutwfc_dict[str(ityp)]:ecutwfc = ecutwfc_dict[str(ityp)]
    if ecutrho < ecutrho_dict[str(ityp)]:ecutrho = ecutrho_dict[str(ityp)]
nq = numpy.zeros(3, numpy.int_)
for ii in range(3):
    norm = numpy.sqrt(numpy.dot(bvec[ii][:], bvec[ii][:]))
    nq[ii] = round(norm / 0.3359385398275)
    print(norm)
nelec = 0.0
for iat in range(nat):
    nelec += valence_dict[str(structure.species[iat])]

print(structure)
#
# scf.in : Charge density
#
with open("scf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'scf\'", file=f)
    print("      outdir = \'./\'", file=f)
    print("  pseudo_dir = \'../pseudo/\'", file=f)
    print("      prefix = \'%s\'" % args[2], file=f)
    print("/", file=f)
    #
    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %f" % ecutwfc, file=f)
    print("     ecutrho = %f" % ecutrho, file=f)
    print(" occupations = \'tetrahedra_opt\'", file=f)
    print("/", file=f)
    #
    print("&ELECTRONS", file=f)
    print("/", file=f)
    #
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii][0],avec[ii][1],avec[ii][2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            structure.species[iat],
            structure.frac_coords[iat][0],structure.frac_coords[iat][1],structure.frac_coords[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)

#
# nscf.in : Dense k grid
#
with open("nscf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'nscf\'", file=f)
    print("      outdir = \'./\'", file=f)
    print("  pseudo_dir = \'../pseudo/\'", file=f)
    print("      prefix = \'%s\'" % args[2], file=f)
    print("/", file=f)
    #
    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %f" % ecutwfc, file=f)
    print("     ecutrho = %f" % ecutrho, file=f)
    print(" occupations = \'tetrahedra_opt\'", file=f)
    print("        nbnd = %d" % int(nelec), file=f)
    print("/", file=f)
    #
    print("&ELECTRONS", file=f)
    print("/", file=f)
    #
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii][0],avec[ii][1],avec[ii][2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            structure.species[iat],
            structure.frac_coords[iat][0],structure.frac_coords[iat][1],structure.frac_coords[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
