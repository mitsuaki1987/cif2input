#!/usr/bin/python3
import os
import pymatgen
from pymatgen.core.periodic_table import get_el_sp
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict
# from sg15 import pseudo_dict, ecutwfc_dict, ecutrho_dict
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
        #
        structure = pymatgen.Structure.from_file(input_file)
        avec = structure.lattice.matrix
        bvec = pymatgen.core.Lattice(structure.lattice.matrix).reciprocal_lattice.matrix
        nat = structure.num_sites
        atom = [str(get_el_sp(iat)) for iat in structure.atomic_numbers]
        typ = set(atom)
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
        spg_analysis = SpacegroupAnalyzer(structure)
        coarse = spg_analysis.get_ir_reciprocal_mesh(mesh=(nk[0], nk[1], nk[2]))
        n_proc = min(4, len(coarse))
        #
        scf_input = "scf_" + prefix + ".in"
        scf_output = "scf_" + prefix + ".out"
        pdos_input = "pdos_" + prefix + ".in"
        pdos_output = "pdos_" + prefix + ".out"
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
                print(" %25.15e %25.15e %25.15e" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
            print("ATOMIC_SPECIES", file=f)
            for ityp in typ:
                print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, pseudo_dict[str(ityp)]), file=f)
            print("ATOMIC_POSITIONS crystal", file=f)
            for iat in range(nat):
                print(" %s %25.15e %25.15e %25.15e" % (
                      atom[iat],
                      structure.frac_coords[iat][0], structure.frac_coords[iat][1], structure.frac_coords[iat][2]),
                      file=f)
            print("K_POINTS automatic", file=f)
            print(" %d %d %d 0 0 0" % (nk[0], nk[1], nk[2]), file=f)
        #
        # Run DFT
        #
        try:
            subprocess.check_call("mpiexec -n %d ~/bin/pw.x -nk %d -in %s > %s"
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
        # proj.in : Read by projwfc.x
        #
        with open(pdos_input, 'w') as f:
            print("&PROJWFC", file=f)
            print(" prefix = \'%s\'" % prefix, file=f)
            print("   emin = %f" % efermi, file=f)
            print("   emax = %f" % efermi, file=f)
            print(" deltae = 0.1", file=f)
            print("/", file=f)
        #
        # Run PDOS
        #
        try:
            subprocess.check_call("mpiexec -n %d ~/bin/projwfc.x -nk %d -in %s > %s"
                                  % (n_proc, n_proc, pdos_input, pdos_output), shell=True)
        except subprocess.CalledProcessError:
            print("PDOS error in ", prefix)
            clean(prefix)
            continue
        #
        clean(prefix)


main()
