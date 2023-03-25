#!/usr/bin/python3
import os
import pymatgen
from pymatgen.core.periodic_table import get_el_sp
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, atomwfc_dict, band_dict
import subprocess
import numpy
import sys


def clean(prefix):
    subprocess.call("rm -rf %s.save %s.dos %s.xml %s.wfc* %s.mix*"
                    % (prefix, prefix, prefix, prefix, prefix), shell=True)


def main():
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
        avec = structure.lattice.matrix
        nat = structure.num_sites
        atom = [str(get_el_sp(iat)) for iat in structure.atomic_numbers]
        typ = sorted(set(atom))
        ntyp = structure.ntypesp
        if nat > 100:
            print("Too many atoms in ", prefix)
            continue
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
        nk = numpy.zeros(3, numpy.int_)
        for ii in range(3):
            norm = numpy.sqrt(numpy.dot(avec[ii][:], avec[ii][:]))
            nk[ii] = round(2.0 * numpy.pi / norm / 0.15)
            if nk[ii] == 0:
                nk[ii] = 1
        #
        scf_input = "scf_" + prefix + ".in"
        scf_output = "scf_" + prefix + ".out"
        nscf_input = "nscf_" + prefix + ".in"
        nscf_output = "nscf_" + prefix + ".out"
        pdos_input = "pdos_" + prefix + ".in"
        pdos_output = "pdos_" + prefix + ".out"
        #
        # Number of valence band
        #
        nbnd = 0
        for iat in atom:
            nbnd += band_dict[iat] + 1
        #
        # SCF file
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
            print("       nspin = 2", file=f)
            print("        nbnd = %d" % nbnd, file=f)
            for ityp in range(ntyp):
                print(" starting_magnetization(%d) = 1.0" % (ityp + 1), file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" diagonalization = \"rmm-davidson\"", file=f)
            print(" mixing_beta = 0.1", file=f)
            print(" conv_thr = %e" % (float(nat)*1.0e-7), file=f)
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
            print(" %d %d %d 0 0 0" % (nk[0], nk[1], nk[2]), file=f)
        #
        # Run DFT (SCF)
        #
        try:
            subprocess.check_call("mpiexec -n %d -of %s ~/bin/pw.x -nk %d -in %s"
                                  % (n_proc, scf_output, n_proc, scf_input), shell=True)
        except subprocess.CalledProcessError:
            print("SCF error in ", prefix)
            clean(prefix)
            continue
        #
        # Unconverged case
        #
        try:
            subprocess.check_call("grep \"convergence has been achieved in\" %s"
                                  % scf_output, shell=True)
        except subprocess.CalledProcessError:
            print("SCF did not converge in ", prefix)
            clean(prefix)
            continue
        #
        # Non-SCF file
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
            print(" occupations = \'tetrahedra_opt\'", file=f)
            print("       nspin = 2", file=f)
            print("        nbnd = %d" % nbnd, file=f)
            for ityp in range(ntyp):
                print(" starting_magnetization(%d) = 1.0" % (ityp + 1), file=f)
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" diagonalization = \"ppcg\"", file=f)
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
            print(" %d %d %d 0 0 0" % (nk[0]*2, nk[1]*2, nk[2]*2), file=f)
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
        # Fermi energy
        #
        efermi = None
        with open(nscf_output, "r") as f:
            output_lines = f.readlines()
        for output_line in output_lines:
            output_words = output_line.split("the Fermi energy is")
            if len(output_words) > 1:
                efermi = float(output_words[1].split()[0])
                #
                # proj.in : Read by projwfc.x
                #
                with open(pdos_input, 'w') as f:
                    print("&PROJWFC", file=f)
                    print(" prefix = \'%s\'" % prefix, file=f)
                    print("   emin = %f" % efermi, file=f)
                    print("   emax = %f" % efermi, file=f)
                    print(" deltae = 0.1", file=f)
                    print("/", file=f)
        if efermi is None:
            print("efermi error in ", prefix)
            clean(prefix)
            continue
        #
        # Run PDOS
        #
        try:
            subprocess.check_call("mpiexec -n %d -of %s ~/bin/projwfc.x -nk %d -in %s"
                                  % (n_proc, pdos_output, n_proc, pdos_input), shell=True)
        except subprocess.CalledProcessError:
            print("PDOS error in ", prefix)
            clean(prefix)
            continue
        #
        # Sum PDOS at each Atom and L
        #
        for ityp in typ:
            nwfc = 1
            for il in atomwfc_dict[ityp]:
                try:
                    subprocess.check_call("mpiexec -n 1 -of %s.pdos_%s%s ~/bin/sumpdos.x %s.pdos_atm*\\(%s\\)_wfc#%d*"
                                          % (prefix, ityp, il[0], prefix, ityp, nwfc), shell=True)
                except subprocess.CalledProcessError:
                    print("Sum-PDOS error in ", prefix)
                    continue
                nwfc += 1
        #
        # Atomwfc dictionary for fermi_proj.x
        #
        pfermi = {ityp: {il: [] for il in atomwfc_dict[ityp][:, 0]} for ityp in typ}
        ii = 0
        for iat in atom:
            for il in atomwfc_dict[iat]:
                for im in range(int(il[1])):
                    ii += 1
                    pfermi[iat][il[0]].append(ii)
        #
        # Fermi surface with atomic projection
        #
        for ityp in typ:
            for il in atomwfc_dict[ityp]:
                with open("fermi_proj.in", 'w') as f:
                    print("&PROJWFC", file=f)
                    print(" prefix = \'%s\'" % prefix, file=f)
                    print("/", file=f)
                    print(len(pfermi[ityp][il]), file=f)
                    for ii in pfermi[ityp][il]:
                        print(" %d" % ii, end="", file=f)
                try:
                    subprocess.check_call("mpiexec -n 1 ~/bin/fermi_proj.x -in fermi_proj.in", shell=True)
                except subprocess.CalledProcessError:
                    print("fermi_proj error in ", prefix)
                    continue
                os.rename("./proj1.frmsf",
                          "./" + prefix + "_" + ityp + il[0] + "_1.frmsf")
                os.rename("./proj2.frmsf",
                          "./" + prefix + "_" + ityp + il[0] + "_2.frmsf")
                # os.rename("./" + prefix + "_proj1.frmsf",
                #           "./" + prefix + "_" + ityp + il[0] + "_1.frmsf")
                # os.rename("./" + prefix + "_proj2.frmsf",
                #           "./" + prefix + "_" + ityp + il[0] + "_2.frmsf")

        clean(prefix)


main()
