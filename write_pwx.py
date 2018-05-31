import os
import pymatgen
from pymatgen.core.periodic_table import get_el_sp


def write_atom(f, avec, typ, nat, pos, atom, pseudo_dict):
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)


def write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel):
    print("  pseudo_dir = \'%s\'" % pseudo_dir, file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
    print("/", file=f)
    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %f" % ecutwfc, file=f)
    print("     ecutrho = %f" % ecutrho, file=f)
    if rel:
        print("     noncolin = .TRUE.", file=f)
        print("     lspinorb = .TRUE.", file=f)
    else:
        print("     noncolin = .FALSE.", file=f)
        print("     lspinorb = .FALSE.", file=f)


def write_pwx(prefix, skp, pseudo_dir, ecutwfc, ecutrho, pseudo_dict, nq, nbnd, rel):
    #
    # Lattice information
    #
    avec = skp["primitive_lattice"]
    pos = skp["primitive_positions"]
    nat = len(skp["primitive_types"])
    atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    typ = set(atom)
    ntyp = len(typ)
    #
    # rx.in : Variation cell optimization
    #
    if not os.path.isfile("rx.in"):
        with open("rx.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'vc-relax\'", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print(" occupations = \'smearing\'", file=f)
            print("    smearing = \'m-p\'", file=f)
            print("     degauss = 0.05", file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" conv_thr = %e" % (float(nat)*1.0e-10), file=f)
            print(" mixing_beta = 0.3", file=f)
            print("/", file=f)
            print("&IONS", file=f)
            print(" ion_dynamics = \"bfgs\"", file=f)
            print("/", file=f)
            print("&CELL", file=f)
            print(" press = 0.0", file=f)
            print(" cell_dynamics = \"bfgs\"", file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # scf.in : Charge density
    #
    if not os.path.isfile("scf.in"):
        with open("scf.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'scf\'", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print(" occupations = \'tetrahedra_opt\'", file=f)
            print("    smearing = \'m-p\'", file=f)
            print("     degauss = 0.05", file=f)
            print("/", file=f)
            #
            print("&ELECTRONS", file=f)
            print(" mixing_beta = 0.3", file=f)
            print(" conv_thr = %e" % (float(nat)*1.0e-10), file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # nscf.in : Dense k grid
    #
    if not os.path.isfile("nscf.in"):
        with open("nscf.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'nscf\'", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print(" occupations = \'tetrahedra_opt\'", file=f)
            print("        nbnd = %d" % nbnd, file=f)
            print("        la2f = .true.", file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
    #
    # band.in : Plot band
    #
    if not os.path.isfile("band.in"):
        with open("band.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'bands\'", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print("        nbnd = %d" % nbnd, file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS crystal", file=f)
            print(len(skp["explicit_kpoints_rel"]), file=f)
            for ik in range(len(skp["explicit_kpoints_rel"])):
                print(" %f %f %f 1.0" % (
                    skp["explicit_kpoints_rel"][ik][0],
                    skp["explicit_kpoints_rel"][ik][1],
                    skp["explicit_kpoints_rel"][ik][2]),
                      file=f)

    #
    # nscf_w.in : Non-scf for wannier90
    #
    if not os.path.isfile("nscf_w.in"):
        with open("nscf_w.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'bands\'", file=f)
            print("  wf_collect = .true.", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print("        nbnd = %d" % nbnd, file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS crystal", file=f)
            print(nq[0]*nq[1]*nq[2], file=f)
            for i0 in range(nq[0]):
                for i1 in range(nq[1]):
                    for i2 in range(nq[2]):
                        print(" %f %f %f %f" % (
                                float(i0)/float(nq[0]),
                                float(i1)/float(nq[1]),
                                float(i2)/float(nq[2]),
                                1.0/float(nq[0]*nq[1]*nq[2])
                        ), file=f)
    #
    # nscf_r.in : Pre-process for respack
    #
    if not os.path.isfile("nscf_r.in"):
        with open("nscf_r.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'nscf\'", file=f)
            print("  wf_collect = .true.", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print(" occupations = \'tetrahedra_opt\'", file=f)
            print("        nbnd = %d" % nbnd, file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nq[0], nq[1], nq[2]), file=f)
    #
    # twin.in : For sctk input
    #
    if not os.path.isfile("twin.in"):
        with open("twin.in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'bands\'", file=f)
            print("  wf_collect = .true.", file=f)
            write_middle(f, pseudo_dir, prefix, nat, ntyp, ecutwfc, ecutrho, rel)
            print("        nbnd = %d" % nbnd, file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print("/", file=f)
            write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
            print("K_POINTS crystal", file=f)
            print(nq[0]*nq[1]*nq[2]*2, file=f)
            #
            # Without Shift
            #
            kvec = [0.0]*3
            for i0 in range(nq[0]):
                kvec[0] = float(i0)/float(nq[0])
                if i0*2 >= nq[0]:
                    kvec[0] += -1.0
                for i1 in range(nq[1]):
                    kvec[1] = float(i1) / float(nq[1])
                    if i1 * 2 >= nq[1]:
                        kvec[1] += -1.0
                    for i2 in range(nq[2]):
                        kvec[2] = float(i2) / float(nq[2])
                        if i2 * 2 >= nq[2]:
                            kvec[2] += -1.0
                        print(" %f %f %f 1.0" % (kvec[0], kvec[1], kvec[2]), file=f)
            #
            # Shifted
            #
            for i0 in range(nq[0]):
                kvec[0] = (float(i0) + 0.5)/float(nq[0])
                if i0 * 2 + 1 >= nq[0]:
                    kvec[0] += -1.0
                for i1 in range(nq[1]):
                    kvec[1] = (float(i1) + 0.5) / float(nq[1])
                    if i1 * 2 + 1 >= nq[1]:
                        kvec[1] += -1.0
                    for i2 in range(nq[2]):
                        kvec[2] = (float(i2) + 0.5) / float(nq[2])
                        if i2 * 2 + 1 >= nq[2]:
                            kvec[2] += -1.0
                        print(" %f %f %f 1.0" % (kvec[0], kvec[1], kvec[2]), file=f)
