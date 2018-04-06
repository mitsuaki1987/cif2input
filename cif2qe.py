#!/usr/bin/python3
import sys
import numpy
import pymatgen
import seekpath
from openmx import omx_pao_dict, omx_pot_dict, omx_radius_dict, omx_valence_dict

args = sys.argv
#
# Pseudo potential type
#
if len(args) > 5:
    if args[5] == "sssp":
        from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict
        print("SSSP id used.")
    else:
        from sg15 import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict
        print("SG15 id used.")
else:
    from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict
    print("SSSP id used.")
#
# Location of pseudo
#
if len(args) > 6:
    pseudo_dir = args[6]
else:
    pseudo_dir = "../pseudo/"
print("Pseudo is at ", pseudo_dir)
#
# CIF parser
#
structure = pymatgen.Structure.from_file(args[1])
structure.remove_oxidation_states()

num2atom = {str(pymatgen.Element(str(spc)).number): str(spc) for spc in structure.species}
#
# Band path and primitive lattice
#
if len(args) > 3:
    reference_distance = float(args[3])
else:
    reference_distance = 0.1
print("  reference_distance : {0}".format(reference_distance))
skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                                    [pymatgen.Element(str(spc)).number for spc in structure.species]),
                                   reference_distance=reference_distance)
#
# Lattice information
#
avec = skp["primitive_lattice"]
bvec = skp["reciprocal_primitive_lattice"]
pos = skp["primitive_positions"]
nat = len(skp["primitive_types"])
atom = [num2atom[str(skp["primitive_types"][iat])] for iat in range(nat)]
typ = set(atom)
ntyp = len(typ)
#
# WFC and Rho cutoff
#
ecutwfc = 0.0
ecutrho = 0.0
for ityp in typ:
    if ecutwfc < ecutwfc_dict[str(ityp)]:
        ecutwfc = ecutwfc_dict[str(ityp)]
    if ecutrho < ecutrho_dict[str(ityp)]:
        ecutrho = ecutrho_dict[str(ityp)]
#
# k and q grid
#
if len(args) > 4:
    reference_distance = float(args[4])
else:
    reference_distance = 0.3359385398275
print("  reference_distance : {0}".format(reference_distance))
nq = numpy.zeros(3, numpy.int_)
for ii in range(3):
    norm = numpy.sqrt(numpy.dot(bvec[ii][:], bvec[ii][:]))
    nq[ii] = round(norm / reference_distance)
    print(norm)
nelec = 0.0
for iat in range(nat):
    nelec += valence_dict[atom[iat]]
#
# Prefix
#
if len(args) > 2:
    prefix = args[2]
else:
    xxx = {}
    for ityp in typ:
        xxx[ityp] = 0
    for iatom in range(nat):
        xxx[atom[iatom]] += 1
    prefix = ""
    for ityp in typ:
        prefix += str(ityp).lower() + str(xxx[ityp])
print("  prefix : {0}".format(prefix))
#
# scf.in : Charge density
#
with open("scf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'scf\'", file=f)
    print("  pseudo_dir = \'%s\'" % pseudo_dir, file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
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
    print(" mixing_beta = 0.3", file=f)
    print("/", file=f)
    #
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
#
# nscf.in : Dense k grid
#
with open("nscf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'nscf\'", file=f)
    print("  pseudo_dir = \'%s\'" % pseudo_dir, file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
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
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
#
# band.in : Plot band
#
with open("band.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'bands\'", file=f)
    print("  pseudo_dir = \'%s\'" % pseudo_dir, file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
    print("/", file=f)
    #
    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %f" % ecutwfc, file=f)
    print("     ecutrho = %f" % ecutrho, file=f)
    print("        nbnd = %d" % int(nelec), file=f)
    print("/", file=f)
    #
    print("&ELECTRONS", file=f)
    print("/", file=f)
    #
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
    print("K_POINTS crystal", file=f)
    print(len(skp["explicit_kpoints_rel"]), file=f)
    for ik in range(len(skp["explicit_kpoints_rel"])):
        print(" %f %f %f 1.0" % (
            skp["explicit_kpoints_rel"][ik][0], skp["explicit_kpoints_rel"][ik][1], skp["explicit_kpoints_rel"][ik][2]),
              file=f)
#
# bands.in : Read by bands.x
#
with open("bands.in", 'w') as f:
    print("&BANDS", file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
    print("       lsym = .false.", file=f)
    print("/", file=f)
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
# proj.in : Read by projwfc.x
#
with open("proj.in", 'w') as f:
    print("&PROJWFC", file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
    print("      emin = ", file=f)
    print("      emax = ", file=f)
    print("    deltae = ", file=f)
    print("/", file=f)
#
# ph.in : Phonon
#
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
# q2r.in : IFC in real space
#
with open("q2r.in", 'w') as f:
    print("&INPUT", file=f)
    print(" fildyn = \'matdyn\'", file=f)
    print("   la2f = .true.", file=f)
    print("   zasr = \'crystal\'", file=f)
    print("  flfrc = \'ifc.dat\'", file=f)
    print("/", file=f)
#
# disp.in : Phonon dispersion
#
with open("disp.in", 'w') as f:
    print("&INPUT", file=f)
    print(" fildyn = \'matdyn\'", file=f)
    print("   la2f = .true.", file=f)
    print("   q_in_cryst_coord = .true.", file=f)
    print("   asr = \'crystal\'", file=f)
    print("  flfrc = \'ifc.dat\'", file=f)
    print("/", file=f)
#
# phdos.in : Phonon DOS
#
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
# nscf_w.in : Non-scf for wannier90
#
with open("nscf_w.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'bands\'", file=f)
    print("  pseudo_dir = \'%s\'" % pseudo_dir, file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
    print("  wf_collect = .true.", file=f)
    print("/", file=f)
    #
    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %f" % ecutwfc, file=f)
    print("     ecutrho = %f" % ecutrho, file=f)
    print("        nbnd = %d" % int(nelec), file=f)
    print("/", file=f)
    #
    print("&ELECTRONS", file=f)
    print("/", file=f)
    #
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
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
# pw2wan.in : PW & wannier90 interface
#
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
# {prefix}.win : wannier90 input
#
with open(prefix + ".win", 'w') as f:
    print("num_bands = %d" % int(nelec), file=f)
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
# pp.in : Plot Kohn-Sham orbitals
#
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
# nscf_r.in : Pre-process for respack
#
with open("nscf_r.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'nscf\'", file=f)
    print("  pseudo_dir = \'%s\'" % pseudo_dir, file=f)
    print("  wf_collect = .true.", file=f)
    print("      prefix = \'%s\'" % prefix, file=f)
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
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0], nq[1], nq[2]), file=f)
#
# respack.in : Input file for
#
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
# rpa.in : Input file for rpa_el.x
#
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
#
# openmx.in : Input file for openmx
#
with open("openmx.in", 'w') as f:
    print("#", file=f)
    print("# File Name", file=f)
    print("#", file=f)
    print("System.CurrrentDirectory    ./", file=f)
    print("System.Name          %s" % prefix, file=f)
    print("level.of.stdout      1", file=f)
    print("level.of.fileout     1", file=f)
    print("data.path            ", file=f)
    print("HS.fileout     off   # on|off, default=off", file=f)
    #
    print("#", file=f)
    print("# Atomic Structure", file=f)
    print("#", file=f)
    print("Species.Number  %d" % ntyp, file=f)
    print("<Definition.of.Atomic.Species", file=f)
    for ityp in typ:
        print(" %s  %s%s-%s  %s  %f" % (
            ityp, ityp, omx_radius_dict[str(ityp)], omx_pao_dict[str(ityp)], omx_pot_dict[str(ityp)],
            pymatgen.Element(ityp).atomic_mass), file=f)
    for ityp in typ:
        print("# proj%s  %s%s-s1p1d1  %s" % (
            ityp, ityp, omx_radius_dict[str(ityp)], omx_pot_dict[str(ityp)]), file=f)
    print("Definition.of.Atomic.Species>", file=f)
    print("Atoms.Number  %d" % nat, file=f)
    print("Atoms.SpeciesAndCoordinates.Unit   FRAC", file=f)
    print("<Atoms.SpeciesAndCoordinates", file=f)
    for iat in range(nat):
        print("%d %s %f %f %f %f %f" % (
            iat+1, atom[iat], pos[iat][0], pos[iat][1], pos[iat][2],
            omx_valence_dict[atom[iat]]*0.5, omx_valence_dict[atom[iat]]*0.5), file=f)
    print("Atoms.SpeciesAndCoordinates>", file=f)
    print("Atoms.UnitVectors.Unit  Ang", file=f)
    print("<Atoms.UnitVectors", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    print("Atoms.UnitVectors>", file=f)
    #
    print("#", file=f)
    print("# SCF or Electronic System", file=f)
    print("#", file=f)
    print("scf.XcType               GGA-PBE", file=f)
    print("scf.SpinPolarization        Off   # On|Off|NC", file=f)
    print("scf.SpinOrbit.Coupling      off   # On|Off", file=f)
    print("<scf.SO.factor", file=f)
    for ityp in typ:
        print(" %s  s 1.0 p 1.0 d 1.0 f 1.0" % ityp, file=f)
    print("scf.SO.factor>", file=f)
    print("scf.maxIter                  40", file=f)
    print("scf.EigenvalueSolver       band        # DC|GDC|Cluster|Band", file=f)
    print("scf.Kgrid              %d %d %d" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    print("scf.Mixing.Type           rmm-diisk    # Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk", file=f)
    print("scf.Init.Mixing.Weight     0.30", file=f)
    print("scf.Min.Mixing.Weight      0.001", file=f)
    print("scf.Max.Mixing.Weight      0.4", file=f)
    print("scf.Mixing.History          30", file=f)
    print("scf.Mixing.StartPulay       6", file=f)
    print("scf.Mixing.EveryPulay       6", file=f)
    print("#scf.Kerker.factor  ", file=f)
    print("scf.criterion             1.0e-6", file=f)
    print("scf.partialCoreCorrection    on", file=f)
    print("scf.restart                off", file=f)
    #
    print("#", file=f)
    print("# DFT+U", file=f)
    print("#", file=f)
    print("scf.Hubbard.U          off", file=f)
    print("scf.Hubbard.Occupation     dual     # onsite|full|dual", file=f)
    print("<Hubbard.U.values", file=f)
    for ityp in typ:
        print(" %s  1s 0.0 1p 0.0 1d 0.0" % ityp, file=f)
    print("Hubbard.U.values>", file=f)
    #
    print("#", file=f)
    print("# 1D FFT", file=f)
    print("#", file=f)
    print("1DFFT.EnergyCutoff       3600", file=f)
    print("1DFFT.NumGridK            900", file=f)
    print("1DFFT.NumGridR            900", file=f)
    #
    print("#", file=f)
    print("# van der Waals", file=f)
    print("#", file=f)
    print("scf.DFTD            off", file=f)
    print("version.dftD         2      # 2|3", file=f)
    print("DFTD.scale6            1   # default=0.75|1.0 (for DFT-D2|DFT-D3)", file=f)
    print("DFTD.scale8       0.7875   # default=0.722|0.7875 (for PBE with zero|bj damping)", file=f)
    print("DFTD.sr6           1.217", file=f)
    print("DFTD.a1           0.4289", file=f)
    print("DFTD.a2           4.4407", file=f)
    print("DFTD.cncut_dftD              40            # default=40 (DFTD.Unit)", file=f)
    print("DFTD.Unit                   Ang", file=f)
    print("DFTD.rcut_dftD             100.0", file=f)
    print("DFTD.d                      20.0", file=f)
    print("DFTD.scale6                 0.75", file=f)
    print("DFTD.IntDirection          1 1 1", file=f)
    print("<DFTD.periodicity", file=f)
    for iat in range(nat):
        print("%d 1" % (iat+1), file=f)
    print("DFTD.periodicity>", file=f)
    #
    print("#", file=f)
    print("# MD and Structure optimization", file=f)
    print("#", file=f)
    print("orbitalOpt.Force.Skip       off", file=f)
    print("MD.Type                     NOMD", file=f)
    print("# NOMD|Opt|NVE|NVT_VS|NVT_VS2|NVT_NH", file=f)
    print("# Opt|DIIS|BFGS|RF|EF|", file=f)
    print("# OptC1|OptC2|OptC3|OptC4|OptC5|RFC5", file=f)
    print("# NEB", file=f)
    print("MD.Opt.DIIS.History     3", file=f)
    print("MD.Opt.StartDIIS      5", file=f)
    print("MD.Opt.EveryDIIS       200", file=f)
    print("MD.maxIter             1", file=f)
    print("MD.Opt.criterion      0.0003", file=f)
    print("MD.Opt.Init.Hessian        Schlegel      # Schlegel|iden", file=f)
    print("<MD.Fixed.XYZ", file=f)
    for iat in range(nat):
        print("%d 0 0 0" % (iat+1), file=f)
    print("MD.Fixed.XYZ>", file=f)
    print("scf.stress.tensor  off  # on|off, default=off", file=f)
    print("<MD.Fixed.Cell.Vectors", file=f)
    print("0 0 0", file=f)
    print("0 0 0", file=f)
    print("0 0 0", file=f)
    print("MD.Fixed.Cell.Vectors>", file=f)
    print("<MD.TempControl", file=f)
    print("3", file=f)
    print("100   2  1000.0  0.0  ", file=f)
    print("400  10   700.0  0.4  ", file=f)
    print("700  40   500.0  0.7  ", file=f)
    print("MD.TempControl>", file=f)
    print("<MD.Init.Velocity", file=f)
    for iat in range(nat):
        print("%d 0.0 0.0 0.0" % (iat+1), file=f)
    print("MD.Init.Velocity>", file=f)
    #
    print("#", file=f)
    print("# Band dispersion", file=f)
    print("#", file=f)
    print("Band.dispersion              off", file=f)
    print("Band.Nkpath  %d" % len(skp["path"]), file=f)
    print("<Band.kpath", file=f)
    for ipath in range(len(skp["path"])):
        start = skp["explicit_segments"][ipath][0]
        final = skp["explicit_segments"][ipath][1] - 1
        print("%d  %f %f %f %f %f %f %s %s" % (
            final - start + 1,
            skp["explicit_kpoints_rel"][start][0],
            skp["explicit_kpoints_rel"][start][1],
            skp["explicit_kpoints_rel"][start][2],
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2],
            skp["explicit_kpoints_labels"][start],
            skp["explicit_kpoints_labels"][final]),
            file=f)
    print("Band.kpath>", file=f)
    #
    print("#", file=f)
    print("# Wannier", file=f)
    print("#", file=f)
    print("Wannier.Func.Calc off", file=f)
    print("Wannier.Func.Num 3", file=f)
    print("Wannier.Outer.Window.Bottom  -1.5", file=f)
    print("Wannier.Outer.Window.Top      7.0", file=f)
    print("Wannier.Inner.Window.Bottom  -1.5", file=f)
    print("Wannier.Inner.Window.Top      1.2", file=f)
    print("Wannier.Initial.Projectors.Unit FRAC", file=f)
    print("<Wannier.Initial.Projectors", file=f)
    print("proj-dxy 0.5 0.5 0.5  0.0 0.0 1.0  1.0 0.0 0.0", file=f)
    print("proj-dxz 0.5 0.5 0.5  0.0 0.0 1.0  1.0 0.0 0.0", file=f)
    print("proj-dyz 0.5 0.5 0.5  0.0 0.0 1.0  1.0 0.0 0.0", file=f)
    print("Wannier.Initial.Projectors>", file=f)
    print("Wannier.Interpolated.Bands             off", file=f)
    print("Wannier.Function.Plot                  off", file=f)
    print("Wannier.Function.Plot.SuperCells      1 1 1", file=f)
    print("Wannier.Minimizing.Scheme             0", file=f)
    print("Wannier.Minimizing.StepLength        2.0", file=f)
    print("Wannier.Minimizing.Secant.Steps       5", file=f)
    print("Wannier.Minimizing.Secant.StepLength 2.0", file=f)
    print("Wannier.Minimizing.Conv.Criterion   1e-8", file=f)
    print("Wannier.Minimizing.Max.Steps         200", file=f)
    print("Wannier.Readin.Overlap.Matrix       on", file=f)
    print("Wannier.Dis.SCF.Max.Steps         200", file=f)
    print("Wannier.Dis.Conv.Criterion        1e-8", file=f)
    print("Wannier.Dis.Mixing.Para           0.5", file=f)
    print("Wannier.MaxShells          12", file=f)
    print("Wannier.Kgrid           8  8  8", file=f)
    #
    print("#", file=f)
    print("# DOS", file=f)
    print("#", file=f)
    print("Dos.fileout      off", file=f)
    print("Dos.Erange       -20.0  20.0", file=f)
    print("Dos.Kgrid     %d %d %d" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
    print("OpticalConductivity.fileout    off", file=f)
    #
    print("#", file=f)
    print("# Energy Decomposition", file=f)
    print("#", file=f)
    print("Energy.Decomposition      off        # on|off", file=f)
    #
    print("#", file=f)
    print("# Orbital optimization", file=f)
    print("#", file=f)
    print("orbitalOpt.Method      off   # Off|Species|Atoms", file=f)
    print("orbitalOpt.Opt.Method     EF     # DIIS|EF", file=f)
    print("orbitalOpt.SD.step       0.001", file=f)
    print("orbitalOpt.HistoryPulay     15", file=f)
    print("orbitalOpt.StartPulay       1", file=f)
    print("orbitalOpt.scf.maxIter     40", file=f)
    print("orbitalOpt.Opt.maxIter    100", file=f)
    print("orbitalOpt.per.MDIter     1000000", file=f)
    print("orbitalOpt.criterion      1.0e-4", file=f)
    print("CntOrb.fileout       off        # on|off", file=f)
    print("Num.CntOrb.Atom    %d" % nat, file=f)
    print("<Atoms.Cont.Orbitals", file=f)
    for iat in range(nat):
        print("%d" % iat, file=f)
    print("Atoms.Cont.Orbitals>", file=f)
    #
    print("#", file=f)
    print("# Order-N", file=f)
    print("#", file=f)
    print("orderN.HoppingRanges    6.0", file=f)
    print("orderN.KrylovH.order    400", file=f)
    print("orderN.Exact.Inverse.S   on    #on|off", file=f)
    print("orderN.KrylovS.order   1600", file=f)
    print("orderN.Recalc.Buffer on #on|off", file=f)
    print("orderN.Expand.Core  on  #on|off", file=f)
    #
    print("#", file=f)
    print("# Electric Field", file=f)
    print("#", file=f)
    print("scf.Electric.Field   0.0 0.0 0.0", file=f)
    #
    print("#", file=f)
    print("# Natural population analysis ", file=f)
    print("#", file=f)
    print("NBO.switch   off # on1|on2", file=f)
    print("NBO.Num.CenterAtoms     5", file=f)
    print("<NBO.CenterAtoms", file=f)
    print("269", file=f)
    print("304", file=f)
    print("323", file=f)
    print("541", file=f)
    print("574", file=f)
    print("NBO.CenterAtoms>", file=f)
    #
    print("#", file=f)
    print("# Magnetic field", file=f)
    print("#", file=f)
    print("scf.Constraint.NC.Spin       off    # on|on2|off", file=f)
    print("scf.Constraint.NC.Spin.v    0.0", file=f)
    print("scf.NC.Zeeman.Spin       off        # on|off", file=f)
    print("scf.NC.Mag.Field.Spin    0.0", file=f)
    print("scf.NC.Zeeman.Orbital      off        # on|off", file=f)
    print("scf.NC.Mag.Field.Orbital   0.0", file=f)
    #
    print("#", file=f)
    print("# ESM", file=f)
    print("#", file=f)
    print("ESM.switch    off    # off, on1=v|v|v, on2=m|v|m, on3=v|v|m, on4=on2+EF", file=f)
    print("ESM.buffer.range      10.0", file=f)
    print("ESM.potential.diff     0.0", file=f)
    print("ESM.wall.position        10.0", file=f)
    print("ESM.wall.height        100.0", file=f)
    #
    print("#", file=f)
    print("# NEB", file=f)
    print("#", file=f)
    print("MD.NEB.Number.Images     10", file=f)
    print("MD.NEB.Spring.Const      0.1", file=f)
    print("<NEB.Atoms.SpeciesAndCoordinates", file=f)
    for iat in range(nat):
        print("%d %s %f %f %f %f %f" % (
            iat+1, atom[iat], pos[iat][0], pos[iat][1], pos[iat][2],
            omx_valence_dict[atom[iat]]*0.5, omx_valence_dict[atom[iat]]*0.5), file=f)
    print("NEB.Atoms.SpeciesAndCoordinates>", file=f)
    #
    print("#", file=f)
    print("# STM Image", file=f)
    print("#", file=f)
    print("partial.charge     off    # on|off", file=f)
    print("partial.charge.energy.window   0.0", file=f)
    #
    print("#", file=f)
    print("# Band Unfolding", file=f)
    print("#", file=f)
    print("Unfolding.Electronic.Band      off    # on|off", file=f)
    print("Unfolding.LowerBound        -10", file=f)
    print("Unfolding.UpperBound          10", file=f)
    print("Unfolding.Nkpoint          4", file=f)
    print("<Unfolding.kpoint", file=f)
    print("K 0.33333333333 0.33333333333 0.0000000000", file=f)
    print("G 0.00000000000 0.00000000000 0.0000000000", file=f)
    print("M 0.50000000000 0.00000000000 0.0000000000", file=f)
    print("K 0.33333333333 0.33333333333 0.0000000000", file=f)
    print("Unfolding.kpoint>", file=f)
    print("Unfolding.desired_totalnkpt    30", file=f)
    print("<Unfolding.ReferenceOrigin", file=f)
    print("0.1 0.2 0.3", file=f)
    print("Unfolding.ReferenceOrigin>", file=f)
    #
    print("#", file=f)
    print("# Draw Kohn-Sham orbital", file=f)
    print("#", file=f)
    print("MO.fileout off", file=f)
    print("num.HOMOs  1", file=f)
    print("num.LUMOs  1", file=f)
    print("MO.Nkpoint  1", file=f)
    print("<MO.kpoint", file=f)
    print("0.0  0.0  0.0", file=f)
    print("MO.kpoint>", file=f)
