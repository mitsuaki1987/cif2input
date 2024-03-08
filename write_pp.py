

def write_pp():
    #
    # bands.in : Read by bands.x
    #
    with open("bands.in", 'w') as f:
        print("&BANDS", file=f)
        print("       lsym = .false.", file=f)
        print("/", file=f)
    #
    # bands.in : Read by bands.x
    #
    with open("plotband.in", 'w') as f:
        print("bands.out", file=f)
        print("-1.0e10 1.0e10", file=f)
        print("band.xmgr", file=f)
    #
    # proj.in : Read by projwfc.x
    #
    with open("proj.in", 'w') as f:
        print("&PROJWFC", file=f)
        print("!      emin = ", file=f)
        print("!      emax = ", file=f)
        print("    deltae = 0.1", file=f)
        print("   filproj = 'bands'", file=f)
        print("/", file=f)
    #
    # dos.in : Read by dos.x
    #
    with open("dos.in", 'w') as f:
        print("&DOS", file=f)
        print("      emin = ", file=f)
        print("      emax = ", file=f)
        print("    deltae = 0.1", file=f)
        print("    bz_sum = \"tetrahedra_opt\"", file=f)
        print("/", file=f)
    #
    # pw2wan.in : PW & wannier90 interface
    #
    with open("pw2wan.in", 'w') as f:
        print("&INPUTPP", file=f)
        print("         outdir = \'./\'", file=f)
        print("         prefix = \'pwscf\'", file=f)
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
    with open("pp.in", 'w') as f:
        print("&INPUTPP ", file=f)
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
    with open("q2r.in", 'w') as f:
        print("&INPUT", file=f)
        print(" fildyn = \'matdyn\'", file=f)
        print("   la2f = .true.", file=f)
        print(" lshift_q = .true.", file=f)
        print("   zasr = \'crystal\'", file=f)
        print("  flfrc = \'ifc.dat\'", file=f)
        print("/", file=f)
