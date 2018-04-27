import os


def write_ph(prefix, nq, ecutwfc, nelec):
    #
    # ph.in : Phonon
    #
    if not os.path.isfile("ph.in"):
        with open("ph.in", 'w') as f:
            print("Phonon", file=f)
            print("&INPUTPH", file=f)
            print("    prefix = \'%s\'" % prefix, file=f)
            print("  lshift_q = .true.", file=f)
            print("     ldisp = .true.", file=f)
            print(" reduce_io = .true.", file=f)
            print("    tr2_ph = 1.0d-15", file=f)
            print(" alpha_mix = 0.3", file=f)
            print("  fildvscf = \'dv\'", file=f)
            print("       nq1 = %d" % nq[0], file=f)
            print("       nq2 = %d" % nq[1], file=f)
            print("       nq3 = %d" % nq[2], file=f)
            print("!  start_q = ", file=f)
            print("!   last_q = ", file=f)
            print("/", file=f)
    #
    # elph.in : Electron-phonon
    #
    if not os.path.isfile("elph.in"):
        with open("elph.in", 'w') as f:
            print("Electron-phonon", file=f)
            print("&INPUTPH", file=f)
            print("          prefix = \'%s\'" % prefix, file=f)
            print("        lshift_q = .true.", file=f)
            print("           ldisp = .true.", file=f)
            print("       reduce_io = .true.", file=f)
            print("             nq1 = %d" % nq[0], file=f)
            print("             nq2 = %d" % nq[1], file=f)
            print("             nq3 = %d" % nq[2], file=f)
            print("!        start_q = ", file=f)
            print("!         last_q = ", file=f)
            print("        fildvscf = \'dv\'", file=f)
            print(" electron_phonon = \'lambda_tetra\'", file=f)
            print("             nk1 = %d" % (nq[0]*4), file=f)
            print("             nk2 = %d" % (nq[1]*4), file=f)
            print("             nk3 = %d" % (nq[2]*4), file=f)
            print("/", file=f)
            print("&INPUTA2F", file=f)
            print(" nfreq = %d" % 100, file=f)
            print("/", file=f)
    #
    # epmat.in : Electron-phonon matrix for SCDFT
    #
    if not os.path.isfile("epmat.in"):
        with open("epmat.in", 'w') as f:
            print("Electron-phonon matrix", file=f)
            print("&INPUTPH", file=f)
            print("          prefix = \'%s\'" % prefix, file=f)
            print("        lshift_q = .true.", file=f)
            print("           ldisp = .true.", file=f)
            print("       reduce_io = .true.", file=f)
            print("             nq1 = %d" % nq[0], file=f)
            print("             nq2 = %d" % nq[1], file=f)
            print("             nq3 = %d" % nq[2], file=f)
            print("!        start_q = ", file=f)
            print("!         last_q = ", file=f)
            print("        fildvscf = \'dv\'", file=f)
            print(" electron_phonon = \'scdft_input\'", file=f)
            print("             nk1 = %d" % nq[0], file=f)
            print("             nk2 = %d" % nq[1], file=f)
            print("             nk3 = %d" % nq[2], file=f)
            print("   elph_nbnd_min = ", file=f)
            print("   elph_nbnd_max = ", file=f)
            print("/", file=f)
    #
    # phdos.in : Phonon DOS
    #
    if not os.path.isfile("phdos.in"):
        with open("phdos.in", 'w') as f:
            print("&INPUT", file=f)
            print(" fildyn = \'matdyn\'", file=f)
            print("   la2f = .true.", file=f)
            print("    dos = .true.", file=f)
            print("    asr = \'crystal\'", file=f)
            print("  flfrc = \'ifc.dat\'", file=f)
            print("    nk1 = %d" % (nq[0]*2), file=f)
            print("    nk2 = %d" % (nq[1]*2), file=f)
            print("    nk3 = %d" % (nq[2]*2), file=f)
            print("   ndos = %d" % 100, file=f)
            print("/", file=f)
    #
    # rpa.in : Input file for rpa_el.x
    #
    if not os.path.isfile("rpa.in"):
        with open("rpa.in", 'w') as f:
            print("&CONTROL", file=f)
            print("      prefix = \'%s\'" % prefix, file=f)
            print("/", file=f)
            print("&SYSTEM", file=f)
            print(" start_q = 1", file=f)
            print("  last_q = 1", file=f)
            print("     nmf = 10", file=f)
            print("  laddxc = .FALSE.", file=f)
            print(" ecutwfc = %f" % ecutwfc, file=f)
            print("     nq1 = %d" % nq[0], file=f)
            print("     nq2 = %d" % nq[1], file=f)
            print("     nq3 = %d" % nq[2], file=f)
            print("/", file=f)
    #
    # scdft.in : Input file for scdft.x
    #
    if not os.path.isfile("scdft.in"):
        with open("scdft.in", 'w') as f:
            print("&CONTROL", file=f)
            print("      prefix = \'%s\'" % prefix, file=f)
            print("/", file=f)
            print("&SYSTEM", file=f)
            print("             temp = 0.1", file=f)
            print("             fbee = 1", file=f)
            print("             lbee = %d" % int(nelec), file=f)
            print("              xic = -1.0", file=f)
            print("              nmf = 10", file=f)
            print("               nx = 100", file=f)
            print("               ne = 50", file=f)
            print("             emin = 1.0e-7", file=f)
            print("             emax = 5.0", file=f)
            print(" electron_maxstep = 100", file=f)
            print("         conv_thr = 1.0e-15", file=f)
            print("/", file=f)
