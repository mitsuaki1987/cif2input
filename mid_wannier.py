#!/usr/bin/python3
import os
import pymatgen
import seekpath
from pymatgen.core.periodic_table import get_el_sp
from sg15 import pseudo_dict, ecutwfc_dict, ecutrho_dict, band_dict, core_dict, wan_dict
import subprocess
import numpy
import sys
import math


def clean(prefix):
    subprocess.call("rm -rf %s.save %s.xml %s.wfc* %s.mix* dir-* bands.out"
                   % (prefix, prefix, prefix, prefix), shell=True)


def main():
    #
    dk_path = 0.1
    dq_grid = 0.27
    # dq_grid = 0.5
    #
    args = sys.argv
    with open(str(args[1]), "r") as f:
        input_list = f.readlines()
    n_proc = int(args[2])
    #
    # Read previous result
    #
    for input_file in input_list:
        #
        input_file = input_file.strip("\n")
        prefix = input_file.split("/")[-1].split(".")[0]
        #
        structure = pymatgen.core.Structure.from_file(input_file)
        #
        #
        #
        frac_coord2 = numpy.array(structure.frac_coords)
        for ipos in range(len(frac_coord2)):
            for iaxis in range(3):
                coord3 = frac_coord2[ipos, iaxis] * 6.0
                if abs(round(coord3) - coord3) < 0.001:
                    frac_coord2[ipos, iaxis] = float(round(coord3)) / 6.0
        #
        skp = seekpath.get_path((structure.lattice.matrix, frac_coord2,
                                 [pymatgen.core.Element(str(spc)).number for spc in structure.species]))
        #
        # Lattice information
        #
        avec = skp["primitive_lattice"]
        bvec = skp["reciprocal_primitive_lattice"]
        atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
        typ = sorted(set(atom))
        nat = len(atom)
        ntyp = len(typ)
        #
        # WFC and Rho cutoff
        #
        ecutwfc = 0.0
        ecutrho = 0.0
        unsupported_element = False
        for ityp in typ:
            if str(ityp) in ecutwfc_dict:
                if ecutwfc < ecutwfc_dict[str(ityp)]:
                    ecutwfc = ecutwfc_dict[str(ityp)]
                if ecutrho < ecutrho_dict[str(ityp)]:
                    ecutrho = ecutrho_dict[str(ityp)]
            else:
                unsupported_element = True
                print("Unsupported element in ", prefix)
                break
        if unsupported_element:
            continue
        #
        # k and q grid
        #  the number of grid proportional to the Height of b
        #  b_i * a_i / |a_i| = 2pi / |a_i| (a_i is perpendicular to other b's)
        #
        nq = numpy.zeros(3, numpy.int_)
        for ii in range(3):
            norm = numpy.sqrt(numpy.dot(avec[ii][:], avec[ii][:]))
            nq[ii] = round(2.0 * numpy.pi / norm / dq_grid)
            if nq[ii] == 0:
                nq[ii] = 1
        #
        scf_input = "scf_" + prefix + ".in"
        scf_output = "scf_" + prefix + ".out"
        band_input = "band_" + prefix + ".in"
        band_output = "band_" + prefix + ".out"
        bands_input = "bands_" + prefix + ".in"
        bands_output = "bands_" + prefix + ".out"
        nscf_input = "nscf_" + prefix + ".in"
        nscf_output = "nscf_" + prefix + ".out"
        respack_input = "respack_" + prefix + ".in"
        respack_output = "respack_" + prefix + ".out"
        #
        # ---------- SCF Calculation --------------------
        #
        with open(scf_input, 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'scf\'", file=f)
            print(" prefix = \'%s\'" % prefix, file=f)
            print("/", file=f)
            print("&SYSTEM", file=f)
            print("       ibrav = 0", file=f)
            print("         nat = %d" % nat, file=f)
            print("        ntyp = %d" % ntyp, file=f)
            print("     ecutwfc = %f" % ecutwfc, file=f)
            print("     ecutrho = %f" % ecutrho, file=f)
            print(" occupations = \'tetrahedra_opt\'", file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" diagonalization = \"cg\"", file=f)
            print(" mixing_beta = 0.3", file=f)
            print(" conv_thr = %e" % (float(nat)*1.0e-7), file=f)
            print(" diagonalization = \"cg\"", file=f)
            print("/", file=f)
            print("CELL_PARAMETERS angstrom", file=f)
            for ii in range(3):
                print(" %25.15e %25.15e %25.15e" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
            print("ATOMIC_SPECIES", file=f)
            for ityp in typ:
                print(" %s %f %s" % (ityp, pymatgen.core.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
            print("ATOMIC_POSITIONS crystal", file=f)
            for iat in range(nat):
                print(" %s %25.15e %25.15e %25.15e" % (
                      atom[iat],
                      structure.frac_coords[iat][0], structure.frac_coords[iat][1], structure.frac_coords[iat][2]),
                      file=f)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
        #
        # Run DFT (SCF)
        #
        try:
            # subprocess.check_call("mpiexec -n %d ~/bin/pw.x -nk %d -in %s > %s"
            #                       % (n_proc, n_proc, scf_input, scf_output), shell=True)
            subprocess.check_call("mpiexec -n %d -of %s ~/bin/pw.x -nk %d -in %s"
                                  % (n_proc, scf_output, n_proc, scf_input), shell=True)
        except subprocess.CalledProcessError:
            print("SCF error in ", prefix)
            clean(prefix)
            continue
        #
        # ----------------- Band Calculation ------------------------
        #
        # Get explicit kpath applicable to bands.x and plotband.x
        #
        kpath = []
        nkpath = []
        #
        for ipath in range(len(skp['path'])):
            dk = numpy.array(skp['point_coords'][skp['path'][ipath][1]]) \
                 - numpy.array(skp['point_coords'][skp['path'][ipath][0]])
            dk = numpy.dot(dk, bvec)
            dknorm = math.sqrt(numpy.dot(dk, dk))
            nkpath0 = max(2, int(dknorm / dk_path))
            nkpath.append(nkpath0)
            for ik in range(nkpath0):
                xkpath = ik / nkpath0
                kpath0 = numpy.array(skp['point_coords'][skp['path'][ipath][0]]) * (1.0 - xkpath) \
                    + numpy.array(skp['point_coords'][skp['path'][ipath][1]]) * xkpath
                kpath.append(kpath0)
            #
            # jump case
            #
            if ipath < len(skp['path']) - 1 and skp['path'][ipath][1] != skp['path'][ipath + 1][0]:
                dk1 = numpy.array(skp['point_coords'][skp['path'][ipath][1]]) - kpath[len(kpath) - 1]
                dk1 = numpy.dot(dk1, bvec)
                dk2 = numpy.array(skp['point_coords'][skp['path'][ipath + 1][0]]) \
                    - numpy.array(skp['point_coords'][skp['path'][ipath][1]])  # jump
                dk2 = numpy.dot(dk2, bvec)
                dk1norm = math.sqrt(numpy.dot(dk1, dk1))
                dk2norm = math.sqrt(numpy.dot(dk2, dk2))
                #
                # If the jump is relatively small, shift the last point closer to the end
                #
                if dk2norm < 10.0 * dk1norm:
                    xkpath = dk2norm / dknorm * 0.09
                    kpath0 = numpy.array(skp['point_coords'][skp['path'][ipath][0]]) * xkpath \
                        + numpy.array(skp['point_coords'][skp['path'][ipath][1]]) * (1.0 - xkpath)
                    kpath[len(kpath) - 1] = kpath0
                kpath.append(numpy.array(skp['point_coords'][skp['path'][ipath][1]]))
        kpath.append(numpy.array(skp['point_coords'][skp['path'][len(skp['path']) - 1][1]]))
        #
        # Number of valence band
        #
        nbnd = 0
        for iat in atom:
            nbnd += band_dict[iat] + 1
        #
        # Band calculation file
        #
        with open(band_input, 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'bands\'", file=f)
            print(" prefix = \'%s\'" % prefix, file=f)
            print("/", file=f)
            print("&SYSTEM", file=f)
            print("       ibrav = 0", file=f)
            print("         nat = %d" % nat, file=f)
            print("        ntyp = %d" % ntyp, file=f)
            print("     ecutwfc = %f" % ecutwfc, file=f)
            print("     ecutrho = %f" % ecutrho, file=f)
            print("        nbnd = %d" % nbnd, file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" diagonalization = \"cg\"", file=f)
            print("/", file=f)
            print("CELL_PARAMETERS angstrom", file=f)
            for ii in range(3):
                print(" %25.15e %25.15e %25.15e" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
            print("ATOMIC_SPECIES", file=f)
            for ityp in typ:
                print(" %s %f %s" % (ityp, pymatgen.core.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
            print("ATOMIC_POSITIONS crystal", file=f)
            for iat in range(nat):
                print(" %s %25.15e %25.15e %25.15e" % (
                      atom[iat],
                      structure.frac_coords[iat][0], structure.frac_coords[iat][1], structure.frac_coords[iat][2]),
                      file=f)
            print("K_POINTS crystal", file=f)
            print(len(kpath), file=f)
            for kpath0 in kpath:
                print(" %f %f %f 1.0" % (kpath0[0], kpath0[1], kpath0[2]), file=f)
        #
        # Run DFT (non-SCF)
        #
        try:
            # subprocess.check_call("mpiexec -n %d ~/bin/pw.x -nk %d -in %s > %s"
            #                       % (n_proc, n_proc, band_input, band_output), shell=True)
            subprocess.check_call("mpiexec -n %d -of %s ~/bin/pw.x -nk %d -in %s"
                                  % (n_proc, band_output, n_proc, band_input), shell=True)
        except subprocess.CalledProcessError:
            print("Band error in ", prefix)
            clean(prefix)
            continue
        #
        # ----------- Bands.x calculation -----------------------------------
        #
        # bands.in : Read by bands.x
        #
        with open(bands_input, 'w') as f:
            print("&BANDS", file=f)
            print(" prefix = \'%s\'" % prefix, file=f)
            print("   lsym = .false.", file=f)
            print("/", file=f)
        #
        # Run bands.x
        #
        try:
            # subprocess.check_call("mpiexec -n %d ~/bin/bands.x -nk %d -in %s > %s"
            #                       % (n_proc, n_proc, bands_input, bands_output), shell=True)
            subprocess.check_call("mpiexec -n %d -of %s ~/bin/bands.x -nk %d -in %s"
                                  % (n_proc, bands_output, n_proc, bands_input), shell=True)
        except subprocess.CalledProcessError:
            print("Bands error in ", prefix)
            clean(prefix)
            continue
        #
        os.rename("./bands.out.gnu", "./" + prefix + ".gnu")
        #
        # ----------- Non-SCF calculation for RESPACK -----------------------
        #
        with open(nscf_input, 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'nscf\'", file=f)
            print(" prefix = \'%s\'" % prefix, file=f)
            print("/", file=f)
            print("&SYSTEM", file=f)
            print("       ibrav = 0", file=f)
            print("         nat = %d" % nat, file=f)
            print("        ntyp = %d" % ntyp, file=f)
            print("     ecutwfc = %f" % ecutwfc, file=f)
            print("     ecutrho = %f" % ecutrho, file=f)
            print("        nbnd = %d" % nbnd, file=f)
            print(" occupations = \'tetrahedra_opt\'", file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" diagonalization = \"cg\"", file=f)
            print("/", file=f)
            print("CELL_PARAMETERS angstrom", file=f)
            for ii in range(3):
                print(" %25.15e %25.15e %25.15e" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
            print("ATOMIC_SPECIES", file=f)
            for ityp in typ:
                print(" %s %f %s" % (ityp, pymatgen.core.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
            print("ATOMIC_POSITIONS crystal", file=f)
            for iat in range(nat):
                print(" %s %25.15e %25.15e %25.15e" % (
                      atom[iat],
                      structure.frac_coords[iat][0], structure.frac_coords[iat][1], structure.frac_coords[iat][2]),
                      file=f)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nq[0], nq[1], nq[2]), file=f)
        #
        # Run DFT (non-SCF)
        #
        try:
            subprocess.check_call("mpiexec -n %d -of %s ~/bin/pw.x -nk %d -in %s"
                                  % (n_proc, nscf_output, n_proc, nscf_input), shell=True)
        except subprocess.CalledProcessError:
            print("Non-SCF error in ", prefix)
            clean(prefix)
            continue
        #
        # --------------- RESPACK preparation --------------------------------------------------
        #
        try:
            subprocess.check_call("python3 ~/bin/qe2respack.py %s.save" % prefix, shell=True)
        except subprocess.CalledProcessError:
            print("QE2RESPACK error in ", prefix)
            clean(prefix)
            continue
        #
        # Number of valence band
        #
        ncore = 0
        for iat in atom:
            ncore += core_dict[iat]
        #
        # Number Wannier
        #
        nwan = 0
        for iat in atom:
            for iwan in wan_dict[iat]:
                if iwan == 's':
                    nwan += 1
                elif iwan == 'p':
                    nwan += 3
                elif iwan == 'd':
                    nwan += 5
        #
        # Energy window
        #
        with open("dir-wfn/dat.eigenvalue", 'r') as f:
            input_list = f.readlines()
        #
        nk0 = int(len(input_list[1:]) / nbnd)
        eig = numpy.array(input_list[1:], dtype=numpy.float_).reshape(nk0, nbnd)
        #
        eig_min = eig[:, ncore:].min()
        eig_max = eig[:, ncore:].max()
        htr2ev = 13.60569228 * 2.0
        eig_min = eig_min * htr2ev - 0.01
        eig_max = eig_max * htr2ev + 0.01
        #
        with open("dir-wfn/dat.bandcalc", 'r') as f:
            input_list = f.readlines()
        ef = float(input_list[1])*htr2ev
        #
        # ---------------- RESPACK Run -------------------------------------------------
        #
        # Input file for RESPACK
        #
        with open(respack_input, 'w') as f:
            print("&PARAM_CHIQW", file=f)
            print("/", file=f)
            print("&PARAM_WANNIER", file=f)
            print("           N_wannier = %d" % nwan, file=f)
            print("     N_initial_guess = %d" % nwan, file=f)
            print(" Lower_energy_window = %f" % eig_min, file=f)
            print(" Upper_energy_window = %f" % eig_max, file=f)
            print("/", file=f)
            for iat in range(nat):
                for iwan in wan_dict[atom[iat]]:
                    if iwan == 's':
                        prjs = ["s"]
                    elif iwan == 'p':
                        prjs = ["px", "py", "pz"]
                    elif iwan == 'd':
                        prjs = ["dxy", "dyz", "dzx", "dx2", "dz2"]
                    else:
                        prjs = []
                    #
                    for prj in prjs:
                        print("%s 0.2 %f %f %f" % (prj,
                                                   structure.frac_coords[iat][0],
                                                   structure.frac_coords[iat][1],
                                                   structure.frac_coords[iat][2]), file=f)
            print("&PARAM_INTERPOLATION", file=f)
            n_sym_points0 = 1
            n_sym_points = []
            for ipath in range(0, len(skp["path"])-1):
                n_sym_points0 += 1
                if skp['path'][ipath][1] != skp['path'][ipath+1][0]:
                    n_sym_points.append(n_sym_points0)
                    n_sym_points0 = 1
            n_sym_points.append(n_sym_points0+1)
            print(" N_sym_points =", end="", file=f)
            for in_sym_point in n_sym_points:
                print(" %d" % in_sym_point, end="", file=f)
            print("", file=f)
            print("/", file=f)
            print("%f %f %f" % (
                skp['point_coords'][skp['path'][0][0]][0],
                skp['point_coords'][skp['path'][0][0]][1],
                skp['point_coords'][skp['path'][0][0]][2]),
                  file=f)
            for ipath in range(len(skp["path"])-1):
                print("%f %f %f" % (
                    skp['point_coords'][skp['path'][ipath][1]][0],
                    skp['point_coords'][skp['path'][ipath][1]][1],
                    skp['point_coords'][skp['path'][ipath][1]][2]),
                      file=f)
                if skp['path'][ipath][1] != skp['path'][ipath+1][0]:
                    print("%f %f %f" % (
                        skp['point_coords'][skp['path'][ipath+1][0]][0],
                        skp['point_coords'][skp['path'][ipath+1][0]][1],
                        skp['point_coords'][skp['path'][ipath+1][0]][2]),
                          file=f)
            print("%f %f %f" % (
                skp['point_coords'][skp['path'][len(skp["path"])-1][1]][0],
                skp['point_coords'][skp['path'][len(skp["path"])-1][1]][1],
                skp['point_coords'][skp['path'][len(skp["path"])-1][1]][2]),
                  file=f)
            print("&PARAM_VISUALIZATION", file=f)
            print("/", file=f)
            print("&PARAM_CALC_INT", file=f)
            print("/", file=f)
        #
        # Run Calc_wannier
        #
        try:
            subprocess.check_call("OMP_NUM_THREADS=48 ~/bin/calc_wannier < %s > %s"
                                  % (respack_input, respack_output), shell=True)
        except subprocess.CalledProcessError:
            print("RESPACK error in ", prefix)
            clean(prefix)
            continue
        #
        # Atomwfc dictionary for fermi_proj.x
        #
        pband = {ityp: {il: [] for il in wan_dict[ityp]} for ityp in typ}
        ii = 0
        for iat in atom:
            for il in wan_dict[iat]:
                if il == 's':
                    n_m = 1
                elif il == 'p':
                    n_m = 3
                elif il == 'd':
                    n_m = 5
                else:
                    n_m = 0
                for im in range(n_m):
                    pband[iat][il].append(ii)
                    ii += 1
        #
        # band.gp : Gnuplot script
        #
        with open(prefix + ".gp", 'w') as f:
            print("EF = %f" % ef, file=f)
            print("Emin = %f" % (eig_min - ef), file=f)
            print("Emax = %f" % (eig_max - ef), file=f)
            print("#", file=f)
            x0 = numpy.linalg.norm(avec[0, :]) * 0.5 / numpy.pi
            print("x1 = 0.0", file=f)
            k0 = 0.0
            for ipath in range(len(skp["path"])):
                dk = numpy.array(skp['point_coords'][skp['path'][ipath][1]]) \
                     - numpy.array(skp['point_coords'][skp['path'][ipath][0]])
                dk = numpy.dot(dk, bvec)
                dknorm = math.sqrt(numpy.dot(dk, dk))
                k0 += dknorm * x0
                print("x%d = %f" % (ipath+2, k0), file=f)
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
            print("set ytics font \'Cmr10,18\'", file=f)
            print("set xtics( \\", file=f)
            label = skp['path'][0][0]
            if label == "GAMMA":
                label = "\\241"
            print("\"%s\" x%d" % (label, 1), end="", file=f)
            for ipath in range(len(skp["path"])-1):
                if skp['path'][ipath][1] == skp['path'][ipath+1][0]:
                    label = skp['path'][ipath][1]
                    if label == "GAMMA":
                        label = "\\241"
                else:
                    label = skp['path'][ipath][1]
                    if label == "GAMMA":
                        label = "\\241"
                    if skp['path'][ipath+1][0] == "GAMMA":
                        label = label + "\\241"
                    else:
                        label = label + skp['path'][ipath+1][0]
                print(", \\\n\"%s\" x%d" % (label, ipath+2), end="", file=f)
            label = skp["path"][len(skp["path"])-1][1]
            if label == "GAMMA":
                label = "\\241"
            print(", \\\n\"%s\" x%d" % (label, len(skp["path"])+1), end="", file=f)
            print(") \\\noffset 0.0, 0.0 font \'Cmr10,18\'", file=f)
            print("#", file=f)
            for ii in range(len(skp["path"])+1):
                print("set arrow from x%d, Emin to x%d, Emax nohead ls 2 front" % (ii+1, ii+1), file=f)
            print("#", file=f)
            print("#", file=f)
            print("set xzeroaxis ls 1", file=f)
            print("#", file=f)
            print("set ylabel \"Energy from {/Cmmi10 E}_F [eV]\" offset - 0.5, 0.0 font \'Cmr10,18\'", file=f)
            print("#", file=f)
            print("set output \"%s.pdf\"" % prefix, file=f)
            print("set terminal pdf color enhanced dashed dl 0.5 size 16.0cm,10.0cm", file=f)
            print("unset key", file=f)
            print("plot[:][Emin:Emax] \\", file=f)
            print("  \"%s.gnu\" u 1:($2-EF) w p, \\" % prefix, file=f)
            print("  \"%s_w.gnu\" u ($1*x%d):($2-EF) w l ls 2" % (prefix, len(skp["path"]) + 1), file=f)
            print("#", file=f)
            print("set output \"%s_fat.pdf\"" % prefix, file=f)
            print("set terminal pdf color enhanced dashed dl 0.5 size 20.0cm, 12.0cm", file=f)
            print("set key outside", file=f)
            print("plot[:][Emin:Emax] \\", file=f)
            for ityp in typ:
                for il in wan_dict[ityp]:
                    print("  \"%s_w.gnu\" u ($1*x%d):($2-EF):((" % (prefix, len(skp["path"]) + 1), end="", file=f)
                    for ii in pband[ityp][il]:
                        print("$%d+" % (ii+3), end="", file=f)
                    print("0)*2.0) every 3 w p ps variable t \"%s %s\", \\" % (ityp, il), file=f)
            print("  \"%s_w.gnu\" u ($1*x%d):($2-EF) w l ls 2 notitle" % (prefix, len(skp["path"]) + 1), file=f)
        #
        # Merge fat band file
        #
        subprocess.check_call("paste dir-wan/dat.iband.fat-* |"
                              + "awk 'NR>2{printf \"%f %f \", $1, $2;"
                              + "for(i=3;i<=NF;i+=3){printf \"%f \", $(i)};print \"\"}' > " + prefix + "_w.gnu",
                              shell=True)
        #
        # Store files
        #
        os.rename("./dir-model/zvo_hr.dat", "./" + prefix + "_hr.dat")
        os.rename("dir-wan/dat.wan-center", "./" + prefix + ".wan-center")

        clean(prefix)


main()
