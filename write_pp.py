import os


def write_pp(prefix):
    #
    # bands.in : Read by bands.x
    #
    if not os.path.isfile("bands.in"):
        with open("bands.in", 'w') as f:
            print("&BANDS", file=f)
            print("      prefix = \'%s\'" % prefix, file=f)
            print("       lsym = .false.", file=f)
            print("/", file=f)
    #
    # proj.in : Read by projwfc.x
    #
    if not os.path.isfile("proj.in"):
        with open("proj.in", 'w') as f:
            print("&PROJWFC", file=f)
            print("      prefix = \'%s\'" % prefix, file=f)
            print("      emin = ", file=f)
            print("      emax = ", file=f)
            print("    deltae = 0.1", file=f)
            print("/", file=f)
    #
    # pw2wan.in : PW & wannier90 interface
    #
    if not os.path.isfile("pw2wan.in"):
        with open("pw2wan.in", 'w') as f:
            print("&INPUTPP", file=f)
            print("         outdir = \'./\'", file=f)
            print("         prefix = \'%s\'" % prefix, file=f)
            print("       seedname = \'%s\'" % prefix, file=f)
            print("      write_mmn = .true.", file=f)
            print("      write_amn = .true.", file=f)
            print("      write_unk = .true.", file=f)
            print("      write_dmn = .true.", file=f)
            print(" spin_component = \'none\'", file=f)
            print("       wan_mode = \'standalone\'", file=f)
            print("/", file=f)
    #
    # pp.in : Plot Kohn-Sham orbitals
    #
    if not os.path.isfile("pp.in"):
        with open("pp.in", 'w') as f:
            print("&INPUTPP ", file=f)
            print("   prefix = \'%s\'" % prefix, file=f)
            print(" plot_num = 7", file=f)
            print("   kpoint = 1", file=f)
            print(" kband(1) = ", file=f)
            print(" kband(2) = ", file=f)
            print("    lsign = .true.", file=f)
            print("/", file=f)
            print("&PLOT  ", file=f)
            print("         iflag = 3", file=f)
            print(" output_format = 5", file=f)
            print("       fileout = \".xsf\"", file=f)
            print("/", file=f)
    #
    # q2r.in : IFC in real space
    #
    if not os.path.isfile("q2r.in"):
        with open("q2r.in", 'w') as f:
            print("&INPUT", file=f)
            print(" fildyn = \'matdyn\'", file=f)
            print("   la2f = .true.", file=f)
            print("   zasr = \'crystal\'", file=f)
            print("  flfrc = \'ifc.dat\'", file=f)
            print("/", file=f)
