import numpy
from pymatgen.core.periodic_table import get_el_sp
from wanorb import wanorb_dict


def write_wannier(skp, nbnd, nq, atomwfc_dict):
    #
    # Lattice information
    #
    avec = skp["primitive_lattice"]
    bvec = skp["reciprocal_primitive_lattice"]
    pos = skp["primitive_positions"]
    nat = len(skp["primitive_types"])
    atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    typ = sorted(set(atom))
    #
    # band.gp : Gnuplot script
    #
    with open("band.gp", 'w') as f:
        print("#set terminal pdf color enhanced \\", file=f)
        print("#dashed dl 0.5 size 8.0cm, 6.0cm", file=f)
        print("#set output \"band.pdf\"", file=f)
        print("#", file=f)
        print("EF = ", file=f)
        print("Emin = ", file=f)
        print("Emax = ", file=f)
        print("#", file=f)
        n_sym_points = 1
        final = 0
        x0 = numpy.linalg.norm(avec[0, :]) * 0.5 / numpy.pi
        print("x%d = %f" % (n_sym_points, x0*skp["explicit_kpoints_linearcoord"][final]), file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start != final:
                n_sym_points += 1
                print("x%d = %f" % (n_sym_points, x0*skp["explicit_kpoints_linearcoord"][start]), file=f)
            n_sym_points += 1
            final = skp["explicit_segments"][ipath][1] - 1
            print("x%d = %f" % (n_sym_points, x0*skp["explicit_kpoints_linearcoord"][final]), file=f)
        print("#", file=f)
        print("set border lw 2", file=f)
        print("#", file=f)
        print("set style line 1 lt 1 lw 2 lc 0 dashtype 2", file=f)
        print("set style line 2 lt 1 lw 2 lc 0", file=f)
        print("set style line 3 lt 1 lw 1 lc 1", file=f)
        print("set style line 4 lt 1 lw 1 lc 2", file=f)
        print("set style line 5 lt 1 lw 1 lc 3", file=f)
        print("set style line 6 lt 1 lw 1 lc 4", file=f)
        print("#", file=f)
        print("set ytics scale 3.0, -0.5 1.0 font \'Cmr10,18\'", file=f)
        print("set xtics( \\", file=f)
        n_sym_points = 1
        final = 0
        label_f = skp["explicit_kpoints_labels"][final]
        if label_f == "GAMMA":
            label_f = "\\241"
        print("\"%s\" x%d" % (label_f, n_sym_points), end="", file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            label_s = skp["explicit_kpoints_labels"][start]
            if label_s == "GAMMA":
                label_s = "\\241"
            label_f = skp["explicit_kpoints_labels"][final]
            if label_f == "GAMMA":
                label_f = "\\241"
            if start != final:
                n_sym_points += 1
                print(", \\\n\"%s%s\" x%d" % (label_f, label_s, n_sym_points), end="", file=f)
            n_sym_points += 1
            final = skp["explicit_segments"][ipath][1] - 1
            label_f = skp["explicit_kpoints_labels"][final]
            if label_f == "GAMMA":
                label_f = "\\241"
            print(", \\\n\"%s\" x%d" % (label_f, n_sym_points), end="", file=f)
        print(") \\\noffset 0.0, 0.0 font \'Cmr10,18\'", file=f)
        print("#", file=f)
        for ii in range(n_sym_points):
            print("set arrow from x%d, Emin to x%d, Emax nohead ls 2 front" % (ii+1, ii+1), file=f)
        print("#", file=f)
        print("set key outside", file=f)
        print("#", file=f)
        print("set xzeroaxis ls 1", file=f)
        print("#", file=f)
        print("set ylabel \"Energy from {/Cmmi10 E}_F [eV]\" offset - 0.5, 0.0 font \'Cmr10,18\'", file=f)
        print("#", file=f)
        n_sym_points = 1
        final = 0
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start == final:
                n_sym_points += 1
                final = skp["explicit_segments"][ipath][1] - 1
            else:
                break
        print("plot[:][Emin:Emax] \\", file=f)
        print("        EF, \\", file=f)
        for ityp in typ:
            for il in range(len(atomwfc_dict[ityp][1])):
                print("        \"" + ityp + atomwfc_dict[ityp][1][il] + ".xmgr\" u 1:2:($3*2) w p ps variable" +
                      " t \"" + ityp + atomwfc_dict[ityp][1][il] + "\", \\",
                      file=f)
        print("        \"wannier_band.dat\" u ($1/%f):($2) w p ls 3, \\" % x0, file=f)
        print("        \"dir-wan/dat.iband\" u ($1*x%d):($2) w l ls 4" % n_sym_points, file=f)
        print("pause -1", file=f)
    #
    # wannier.win : wannier90 input
    #
    with open("wannier.win", 'w') as f:
        print("num_bands = %d" % nbnd, file=f)
        print(" num_wann = ", file=f)
        print("", file=f)
        print(" dis_win_min = ", file=f)
        print(" dis_win_max = ", file=f)
        print("dis_froz_min = ", file=f)
        print("dis_froz_max = ", file=f)
        print("", file=f)
        print("begin projections", file=f)
        print("end projections", file=f)
        print("!site_symmetry = .true.", file=f)
        print("", file=f)
        print("write_hr = .true.", file=f)
        print("bands_plot = .true.", file=f)
        print("wannier_plot = .true.", file=f)
        print("", file=f)
        print("wannier_plot_supercell = 3", file=f)
        print("begin kpoint_path", file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            final = skp["explicit_segments"][ipath][1] - 1
            print("%s %f %f %f %s %f %f %f" % (
                skp["explicit_kpoints_labels"][start],
                skp["explicit_kpoints_rel"][start][0],
                skp["explicit_kpoints_rel"][start][1],
                skp["explicit_kpoints_rel"][start][2],
                skp["explicit_kpoints_labels"][final],
                skp["explicit_kpoints_rel"][final][0],
                skp["explicit_kpoints_rel"][final][1],
                skp["explicit_kpoints_rel"][final][2]),
                  file=f)
        print("end kpoint_path", file=f)
        print("", file=f)
        print("mp_grid = %d %d %d" % (nq[0], nq[1], nq[2]), file=f)
        print("", file=f)
        print("begin unit_cell_cart", file=f)
        print("Ang", file=f)
        for ii in range(3):
            print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
        print("end unit_cell_cart", file=f)
        print("", file=f)
        print("begin atoms_frac", file=f)
        for iat in range(nat):
            print(" %s %f %f %f" % (
                atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
        print("end atoms_frac", file=f)
        print("", file=f)
        print("begin kpoints", file=f)
        for i0 in range(nq[0]):
            for i1 in range(nq[1]):
                for i2 in range(nq[2]):
                    print(" %f %f %f" % (
                            float(i0)/float(nq[0]),
                            float(i1)/float(nq[1]),
                            float(i2)/float(nq[2])
                    ), file=f)
        print("end kpoints", file=f)
    #
    # respack.in : Input file for RESPACK
    #
    with open("respack.in", 'w') as f:
        print("&PARAM_CHIQW", file=f)
        print("          Num_freq_grid = 1", file=f)
        print("!          Ecut_for_eps = ", file=f)
        print("               flg_cRPA = 1", file=f)
        print("! MPI_num_proc_per_qcomm = 1", file=f)
        print("!          MPI_num_qcomm = 1", file=f)
        print("!          flg_calc_type = 2", file=f)
        print("!               n_calc_q = 1", file=f)
        print("/", file=f)
        print("&PARAM_WANNIER", file=f)
        print("           N_wannier = ", file=f)
        print("     N_initial_guess = ", file=f)
        print(" Lower_energy_window = ", file=f)
        print(" Upper_energy_window = ", file=f)
        print("flg_initial_guess_direc = 1", file=f)
        print("!   set_inner_window =.true.", file=f)
        print("! Lower_inner_window = ", file=f)
        print("! Upper_inner_window = ", file=f)
        print("/", file=f)
        for iat in range(nat):
            for prj in wanorb_dict[atom[iat]]:
                print("%s 0.2 %f %f %f  1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0 !%s%d" %
                      (prj, pos[iat][0], pos[iat][1], pos[iat][2],
                       atom[iat], iat+1), file=f)
        print("&PARAM_INTERPOLATION", file=f)
        n_sym_points = 1
        final = 0
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start == final:
                n_sym_points += 1
                final = skp["explicit_segments"][ipath][1] - 1
            else:
                break
        print(" N_sym_points = %d" % n_sym_points, file=f)
        print("!       dense = %d, %d, %d" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
        print("/", file=f)
        final = 0
        print("%f %f %f" % (
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2]),
              file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start == final:
                final = skp["explicit_segments"][ipath][1] - 1
                print("%f %f %f" % (
                    skp["explicit_kpoints_rel"][final][0],
                    skp["explicit_kpoints_rel"][final][1],
                    skp["explicit_kpoints_rel"][final][2]),
                    file=f)
            else:
                break
        print("&PARAM_VISUALIZATION", file=f)
        print("! flg_vis_wannier = 1,", file=f)
        print("       ix_vis_min = -1,", file=f)
        print("       ix_vis_max = 2,", file=f)
        print("       iy_vis_min = -1,", file=f)
        print("       iy_vis_max = 2,", file=f)
        print("       iz_vis_min = -1,", file=f)
        print("       iz_vis_max = 2", file=f)
        print("/", file=f)
        print("&PARAM_CALC_INT", file=f)
        print("  calc_ifreq = 1", file=f)
        print(" ix_intJ_min = 0", file=f)
        print(" ix_intJ_max = 0", file=f)
        print(" iy_intJ_min = 0", file=f)
        print(" iy_intJ_max = 0", file=f)
        print(" iz_intJ_min = 0", file=f)
        print(" iz_intJ_max = 0", file=f)
        print("/", file=f)
    #
    # disp.in : Phonon dispersion
    #
    with open("disp.in", 'w') as f:
        print("&INPUT", file=f)
        print(" flfrc = \'ifc.dat\'", file=f)
        print(" fldos = \' \'", file=f)
        print(" flfrq = \'matdyn.freq\'", file=f)
        print(" flvec = \'matdyn.modes\'", file=f)
        print(" fleig = \' \'", file=f)
        print(" fldyn = \' \'", file=f)
        print(" fltau = \' \'", file=f)
        print("   la2f = .true.", file=f)
        print("   q_in_cryst_coord = .true.", file=f)
        print("   asr = \'crystal\'", file=f)
        print("/", file=f)
        print(len(skp["explicit_kpoints_rel"]), file=f)
        for ik in range(len(skp["explicit_kpoints_rel"])):
            print(" %f %f %f 1.0" % (
                skp["explicit_kpoints_rel"][ik][0],
                skp["explicit_kpoints_rel"][ik][1],
                skp["explicit_kpoints_rel"][ik][2]),
                  file=f)
    #
    # respack.in : Input file for RESPACK
    #
    with open("dcore.ini", 'w') as f:
        print("[model]", file=f)
        print("lattice = wannier90", file=f)
        print("ncor = ", file=f)
        print("nelec = ", file=f)
        print("norb = ", file=f)
        print("seedname = wannier", file=f)
        print("equiv = None", file=f)
        print("bvec = [(%f, %f, %f)," % (bvec[0][0], bvec[0][1], bvec[0][2]), file=f)
        print("        (%f, %f, %f)," % (bvec[1][0], bvec[1][1], bvec[1][2]), file=f)
        print("        (%f, %f, %f)]" % (bvec[2][0], bvec[2][1], bvec[2][2]), file=f)
        print("spin_orbit = False", file=f)
        print("interaction = respack", file=f)
        print("density_density = False", file=f)
        print("kanamori = None", file=f)
        print("slater_f = None", file=f)
        print("slater_uj = None", file=f)
        print("non_colinear = False", file=f)
        print("", file=f)
        print("[system]", file=f)
        print("beta = 40.0", file=f)
        print("n_iw = 2048", file=f)
        print("n_tau = 10000", file=f)
        print("fix_mu = False", file=f)
        print("mu = 0.0", file=f)
        print("nk0 = %d" % (nq[0]*4), file=f)
        print("nk1 = %d" % (nq[1]*4), file=f)
        print("nk2 = %d" % (nq[2]*4), file=f)
        print("prec_mu = 0.0001", file=f)
        print("with_dc = True", file=f)
        print("perform_tail_fit = False", file=f)
        print("fit_max_moment = 2", file=f)
        print("fit_min_w = 5.0", file=f)
        print("fit_max_w = 10.0", file=f)
        print("n_l = 0", file=f)
        print("", file=f)
        print("[impurity_solver]", file=f)
        print("verbosity{int} = 10", file=f)
        print("#name = TRIQS/hubbard-I", file=f)
        print("#name = TRIQS/cthyb", file=f)
        print("#n_cycles{int} = 5000", file=f)
        print("#n_warmup_cycles{int} = 5000", file=f)
        print("#length_cycle{int} = 50", file=f)
        print("name = ALPS/cthyb", file=f)
        print("thermalization_time{int} = 60", file=f)
        print("max_time{int} = 120", file=f)
        print("", file=f)
        print("[control]", file=f)
        print("max_step = 100", file=f)
        print("sigma_mix = 0.5", file=f)
        print("restart = False", file=f)
        print("", file=f)
        print("[tool]", file=f)
        print("nk_line = 20", file=f)
        final = 0
        n_sym_points = 1
        print("knode = [(%s, %f, %f, %f)" % (
            skp["explicit_kpoints_labels"][final],
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2]),
              file=f, end="")
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start == final:
                n_sym_points += 1
                final = skp["explicit_segments"][ipath][1] - 1
                print(",\n         (%s, %f, %f, %f)" % (
                    skp["explicit_kpoints_labels"][final],
                    skp["explicit_kpoints_rel"][final][0],
                    skp["explicit_kpoints_rel"][final][1],
                    skp["explicit_kpoints_rel"][final][2]),
                    file=f, end="")
            else:
                break
        print("]", file=f)
        print("nnode = %d" % n_sym_points, file=f)
        print("omega_min = -1", file=f)
        print("omega_max = 1", file=f)
        print("Nomega = 100", file=f)
        print("broadening = 0.1", file=f)
        print("eta = 0.0", file=f)
        print("omega_pade = 5.0", file=f)
        print("omega_check = 5.0", file=f)
