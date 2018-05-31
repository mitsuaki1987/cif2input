import os
import numpy
from pymatgen.core.periodic_table import get_el_sp


def write_wannier(prefix, skp, nbnd, nq):
    #
    # Lattice information
    #
    avec = skp["primitive_lattice"]
    pos = skp["primitive_positions"]
    nat = len(skp["primitive_types"])
    atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    #
    # band.gp : Gnuplot script
    #
    if not os.path.isfile("band.gp"):
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
            print("unset key", file=f)
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
            print("        \"bands.out.gnu\" u 1:($2-EF) w p ls 3, \\", file=f)
            print("        \"%s_band.dat\" u ($1/%f):($2-EF) w p ls 3, \\" % (prefix, x0), file=f)
            print("        \"dir-wan/dat.iband\" u ($1*x%d):($2-EF) w l ls 4" % n_sym_points, file=f)
            print("pause -1", file=f)
    #
    # {prefix}.win : wannier90 input
    #
    if not os.path.isfile(prefix + ".win"):
        with open(prefix + ".win", 'w') as f:
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
    if not os.path.isfile("respack.in"):
        with open("respack.in", 'w') as f:
            print("&PARAM_CHIQW", file=f)
            print("          Num_freq_grid = 1", file=f)
            print("!          Ecut_for_eps = ", file=f)
            print("               flg_cRPA = 1", file=f)
            print(" MPI_num_proc_per_qcomm = 1", file=f)
            print("          MPI_num_qcomm = 1", file=f)
            print("          flg_calc_type = 2", file=f)
            print("               n_calc_q = 1", file=f)
            print("/", file=f)
            print("1", file=f)
            print("&PARAM_WANNIER", file=f)
            print("           N_wannier = ", file=f)
            print("     N_initial_guess = ", file=f)
            print(" Lower_energy_window = ", file=f)
            print(" Upper_energy_window = ", file=f)
            print("!   set_inner_window =.true.", file=f)
            print("! Lower_inner_window = ", file=f)
            print("! Upper_inner_window = ", file=f)
            print("/", file=f)
            print("", file=f)
            n_sym_points = 1
            final = 0
            for ipath in range(len(skp["path"])):
                start = skp["explicit_segments"][ipath][0]
                if start == final:
                    n_sym_points += 1
                    final = skp["explicit_segments"][ipath][1] - 1
                else:
                    break
            print("&PARAM_INTERPOLATION", file=f)
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
    if not os.path.isfile("disp.in"):
        with open("disp.in", 'w') as f:
            print("&INPUT", file=f)
            print(" fildyn = \'matdyn\'", file=f)
            print("   la2f = .true.", file=f)
            print("   q_in_cryst_coord = .true.", file=f)
            print("   asr = \'crystal\'", file=f)
            print("  flfrc = \'ifc.dat\'", file=f)
            print("/", file=f)
            print(len(skp["explicit_kpoints_rel"]), file=f)
            for ik in range(len(skp["explicit_kpoints_rel"])):
                print(" %f %f %f 1.0" % (
                    skp["explicit_kpoints_rel"][ik][0],
                    skp["explicit_kpoints_rel"][ik][1],
                    skp["explicit_kpoints_rel"][ik][2]),
                      file=f)
