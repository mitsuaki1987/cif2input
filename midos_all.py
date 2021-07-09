#!/usr/bin/python3
import os
import pymatgen
import seekpath
from pymatgen.core.periodic_table import get_el_sp
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict
from xml.etree import ElementTree
import subprocess
import numpy
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys


def clean(prefix):
    subprocess.call("rm -rf %s.save %s.dos %s.xml %s.wfc* %s.mix*"
                    % (prefix, prefix, prefix, prefix, prefix), shell=True)


def main():
    #
    args = sys.argv
    with open(str(args[1]), "r") as f:
        input_list = f.readlines()
    #
    # Read previous result
    #
    for input_file in input_list:
        #
        input_file = input_file.strip("\n")
        prefix = input_file.split("/")[-1].split(".")[0]
        dos_file = 0
        if os.path.isfile("dos_" + prefix + ".dat"):
            dos_file = os.path.getsize("dos_" + prefix + ".dat")
        if dos_file != 0:
            print("already done", prefix)
            continue
        #
        structure = pymatgen.core.Structure.from_file(input_file)
        structure.remove_oxidation_states()
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
        pos = skp["primitive_positions"]
        nat = len(skp["primitive_types"])
        atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
        typ = set(atom)
        ntyp = len(typ)
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
        # k grid
        #
        nk = numpy.zeros(3, numpy.int_)
        for ii in range(3):
            norm = numpy.sqrt(numpy.dot(bvec[ii][:], bvec[ii][:]))
            nk[ii] = round(norm / 0.16)
            if nk[ii] == 0:
                nk[ii] = 1
        #
        # Number of k in IBZ
        #
        structure2 = pymatgen.core.Structure(skp["primitive_lattice"],
                                             skp["primitive_types"],
                                             skp["primitive_positions"])
        spg_analysis = SpacegroupAnalyzer(structure2)
        coarse = spg_analysis.get_ir_reciprocal_mesh(mesh=(nk[0], nk[1], nk[2]))
        n_proc = min(28, len(coarse))
        #
        scf_input = "scf_" + prefix + ".in"
        scf_output = "scf_" + prefix + ".out"
        dos_input = "dos_" + prefix + ".in"
        dos_output = "dos_" + prefix + ".out"
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
            print("/", file=f)
            print("&ELECTRONS", file=f)
            print(" mixing_beta = 0.3", file=f)
            print(" conv_thr = %e" % (float(nat)*1.0e-7), file=f)
            print("/", file=f)
            print("CELL_PARAMETERS angstrom", file=f)
            for ii in range(3):
                print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
            print("ATOMIC_SPECIES", file=f)
            for ityp in typ:
                print(" %s %f %s" % (ityp, pymatgen.core.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
            print("ATOMIC_POSITIONS crystal", file=f)
            for iat in range(nat):
                print(" %s %f %f %f" % (
                    atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nk[0], nk[1], nk[2]), file=f)
        #
        # Run DFT
        #
        try:
            subprocess.check_call("mpirun -hostfile $PBS_NODEFILE -np %d ~/bin/pw.x -nk %d -in %s > %s"
                                  % (n_proc, n_proc, scf_input, scf_output), shell=True)
        except subprocess.CalledProcessError:
            print("SCF error in ", prefix)
            clean(prefix)
            continue
        #
        # Extract DOS
        #
        xmlfile = os.path.join(prefix + ".save/", 'data-file-schema.xml')
        tree = ElementTree.parse(xmlfile)
        root = tree.getroot()
        child = root.find('output').find('band_structure')
        efermi = float(child.find('fermi_energy').text) * 13.60569228 * 2.0
        #
        # DOS file
        #
        with open(dos_input, 'w') as f:
            print("&DOS", file=f)
            print("    prefix = \'%s\'" % prefix, file=f)
            print("      emin = %f" % efermi, file=f)
            print("      emax = %f" % efermi, file=f)
            print("    deltae = 0.1", file=f)
            print("    bz_sum = \"tetrahedra_opt\"", file=f)
            print("/", file=f)
        #
        # Run DOS
        #
        try:
            subprocess.check_call("mpirun -np 1 ~/bin/dos.x -in %s > %s" % (dos_input, dos_output), shell=True)
        except subprocess.CalledProcessError:
            print("DOS error in ", prefix)
            clean(prefix)
            continue
        #
        with open(prefix + ".dos", "r") as f:
            f.readline()
            line = f.readline()
            dos = float(line.split()[1]) / float(nat)
        #
        with open("dos_" + prefix + ".dat", "w") as f:
            print(dos, file=f)
        #
        clean(prefix)


main()
