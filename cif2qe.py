#!/usr/bin/python3
import sys
import numpy
import pymatgen
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict
import seekpath

args = sys.argv

structure = pymatgen.Structure.from_file(args[1])
structure.remove_oxidation_states()

num2atom = {str(pymatgen.Element(str(spc)).number):str(spc) for spc in structure.species}

if len(args) > 3:
    reference_distance = float(args[3])
else:
    reference_distance = 0.1
print("  reference_distance : {0}".format(reference_distance))
skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                         [pymatgen.Element(str(spc)).number for spc in structure.species]),
                                   reference_distance=reference_distance)

avec = skp["primitive_lattice"]
bvec = skp["reciprocal_primitive_lattice"]
pos = skp["primitive_positions"]
nat = len(skp["primitive_types"])
atom = [num2atom[str(skp["primitive_types"][iat])] for iat in range(nat)]
typ = set(atom)
ntyp = len(typ)
#
ecutwfc = 0.0
ecutrho = 0.0
for ityp in typ:
    if ecutwfc < ecutwfc_dict[str(ityp)]:ecutwfc = ecutwfc_dict[str(ityp)]
    if ecutrho < ecutrho_dict[str(ityp)]:ecutrho = ecutrho_dict[str(ityp)]
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
for iat in range(nat): nelec += valence_dict[atom[iat]]
#
if len(args) > 2:
    prefix = args[2]
else:
    xxx = {}
    for ityp in typ: xxx[ityp] = 0
    for iatom in range(nat): xxx[atom[iatom]] += 1
    prefix = ""
    for ityp in typ:
        prefix += str(ityp) + str(xxx[ityp])
print("  prefix : {0}".format(prefix))
#
# scf.in : Charge density
#
with open("scf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'scf\'", file=f)
    print("  pseudo_dir = \'../pseudo/\'", file=f)
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
            atom[iat], pos[iat][0],pos[iat][1],pos[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
#
# nscf.in : Dense k grid
#
with open("nscf.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'nscf\'", file=f)
    print("  pseudo_dir = \'../pseudo/\'", file=f)
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
        print(" %f %f %f" % (avec[ii][0],avec[ii][1],avec[ii][2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0],pos[iat][1],pos[iat][2]), file=f)
    #
    print("K_POINTS automatic", file=f)
    print(" %d %d %d 0 0 0" % (nq[0]*4, nq[1]*4, nq[2]*4), file=f)
#
# band.in : Plot band
#
with open("band.in", 'w') as f:
    print("&CONTROL", file=f)
    print(" calculation = \'bands\'", file=f)
    print("  pseudo_dir = \'../pseudo/\'", file=f)
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
        print(" %f %f %f" % (avec[ii][0],avec[ii][1],avec[ii][2]), file=f)
    #
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
    #
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0],pos[iat][1],pos[iat][2]), file=f)
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
    print("!       lsym = .true.", file=f)
    print("/", file=f)
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
            skp["explicit_kpoints_rel"][start][0], skp["explicit_kpoints_rel"][start][1],  skp["explicit_kpoints_rel"][start][2],
            skp["explicit_kpoints_labels"][final],
            skp["explicit_kpoints_rel"][final][0], skp["explicit_kpoints_rel"][final][1],  skp["explicit_kpoints_rel"][final][2]),
              file=f)
    print("end kpoint_path", file=f)
    print("", file=f)
    print("mp_grid = 4 4 4", file=f)
    print("", file=f)
    print("begin unit_cell_cart", file=f)
    print("Ang", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii][0],avec[ii][1],avec[ii][2]), file=f)
    print("end unit_cell_cart", file=f)
    print("", file=f)
    print("begin atoms_frac", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0],pos[iat][1],pos[iat][2]), file=f)
    print("end atoms_frac", file=f)
    print("", file=f)
    print("begin kpoints", file=f)
    print("end kpoints", file=f)
