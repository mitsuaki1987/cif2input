import pymatgen
from pymatgen.core.periodic_table import get_el_sp


def write_atom(f, avec, typ, nat, pos, atom, pseudo_dict):
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %.12f %.12f %.12f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.core.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %.12f %.12f %.12f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)


def write_head(f, calculation, nat, ntyp, ecutwfc, ecutrho, rel):
    print("&CONTROL", file=f)
    print(" calculation = \'%s\'" % calculation, file=f)
    print("     disk_io = \'low\'", file=f)
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


def write_pwx(avec, atom, pos, ecutwfc, ecutrho, pseudo_dict, nq, nbnd, rel, kpath):
    #
    # Lattice information
    #
    nat = len(atom)
    typ = sorted(set(atom))
    ntyp = len(typ)
    #
    # rx.in : Variation cell optimization
    #
    with open("rx.in", 'w') as f:
        write_head(f, "vc-relax", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" conv_thr = %e" % (float(nat)*1.0e-8), file=f)
        print(" mixing_beta = 0.1", file=f)
        print(" scf_must_converge = .false.", file=f)
        print(" electron_maxstep = 30", file=f)
        print(" diagonalization = \"cg\"", file=f)
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
    with open("scf.in", 'w') as f:
        write_head(f, "scf", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("/", file=f)
        #
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"rmm-davidson\"", file=f)
        print(" mixing_beta = 0.1", file=f)
        print(" conv_thr = %e" % (float(nat)*1.0e-10), file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # nscf.in : Dense k grid
    #
    with open("nscf.in", 'w') as f:
        write_head(f, "nscf", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("        nbnd = %d" % nbnd, file=f)
        print("        la2f = .true.", file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
    #
    # nscf_p.in : Phonon pre-process
    #
    with open("nscf_p.in", 'w') as f:
        write_head(f, "nscf", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # nscf_pd.in : Electron-Phonon pre-process
    #
    with open("nscf_pd.in", 'w') as f:
        write_head(f, "nscf", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
    #
    # nscf_pc.in : Electron-Phonon matrix pre-process
    #
    with open("nscf_pc.in", 'w') as f:
        write_head(f, "nscf", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nq[0], nq[1], nq[2]), file=f)
    #
    # band.in : Plot band
    #
    with open("band.in", 'w') as f:
        write_head(f, "bands", nat, ntyp, ecutwfc, ecutrho, rel)
        print("        nbnd = %d" % nbnd, file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS crystal", file=f)
        print(len(kpath), file=f)
        for kpath0 in kpath:
            print(" %f %f %f 1.0" % (kpath0[0], kpath0[1], kpath0[2]), file=f)
    #
    # nscf_w.in : Non-scf for wannier90
    #
    with open("nscf_w.in", 'w') as f:
        write_head(f, "bands", nat, ntyp, ecutwfc, ecutrho, rel)
        print("        nbnd = %d" % nbnd, file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
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
    with open("nscf_r.in", 'w') as f:
        write_head(f, "nscf", nat, ntyp, ecutwfc, ecutrho, rel)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("        nbnd = %d" % nbnd, file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom, pseudo_dict)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nq[0], nq[1], nq[2]), file=f)
    #
    # twin.in : For sctk input
    #
    with open("twin.in", 'w') as f:
        write_head(f, "bands", nat, ntyp, ecutwfc, ecutrho, rel)
        print("        nbnd = %d" % nbnd, file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print(" diagonalization = \"ppcg\"", file=f)
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
