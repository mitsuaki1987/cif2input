#!/usr/bin/python3
import os
import pymatgen
import seekpath
from pymatgen.core.periodic_table import get_el_sp
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import physbo
import subprocess
import numpy


def load_descriptor():
    with open("desc.dat", "r") as f:
        lines = f.readlines()
        filename = []
        descriptor = []
        for line in lines:
            filename.append(str(line.split()[0]))
            descriptor.append(numpy.array(line.split()[1:], dtype=numpy.float_))

    return numpy.array(descriptor), filename


def load_result(num_action):
    action = []
    result = []
    for i_action in range(num_action):
        if os.path.isfile(str(num_action) + "/dos.dat"):
            action.append(int(i_action))
            with open(str(num_action) + "/dos.dat", 'r') as f:
                result.append(float(f.readline()))

    return numpy.array(action), numpy.array(result)


def qsub_action(file_name, i_action):
    structure = pymatgen.core.Structure.from_file(file_name)
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
    # job file
    #
    with open("scf.sh", 'w') as f:
        print("#!/bin/sh -e", file=f)
        print("#PBS -l nodes=1:ppn=28", file=f)
        print("#PBS -l walltime=8:00:00", file=f)
        print("#PBS -n", file=f)
        print("source ~/.bashrc", file=f)
        print("mkdir $PBS_O_WORKDIR/%d" % i_action, file=f)
        print("cd $PBS_O_WORKDIR/%d" % i_action, file=f)
        #
        # SCF input file
        #
        print("cat > scf.in << EOF", file=f)
        print("&CONTROL", file=f)
        print(" calculation = \'scf\'", file=f)
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
        print("EOF", file=f)
        #
        print("mpijob -n %d ~/bin/pw.x -nk %d -in scf.in > scf.out"
              % (n_proc, n_proc), file=f)
        print("ef=`grep Fermi scf.out| awk '{print $5}'`", file=f)
        #
        # DOS input file
        #
        print("cat > dos.in << EOF", file=f)
        print("&DOS", file=f)
        print("      emin = ${ef}", file=f)
        print("      emax = ${ef}", file=f)
        print("    deltae = 0.1", file=f)
        print("    bz_sum = \"tetrahedra_opt\"", file=f)
        print("/", file=f)
        print("EOF", file=f)
        #
        print("mpijob -n %d ~/bin/dos.x -in dos.in > dos.out"
              % n_proc, file=f)
        #
        # DOS per atom
        #
        print("awk \'NR==2{print $2/%d}\' pwscf.dos > dos.dat" % nat, file=f)
    #
    # Submit batch job
    #
    subprocess.call("qsub scf.sh", shell=True)


def main():
    descriptor, filename = load_descriptor()
    descriptor = physbo.misc.centering(descriptor)
    policy = physbo.search.discrete.policy(test_X=descriptor)
    policy.set_seed(1)
    #
    # Read previous result
    #
    action, result = load_result(len(descriptor))
    if len(action) == 0:
        action = policy.random_search(max_num_probes=1, num_search_each_probe=5, simulator=None)
    else:
        policy.write(action, result)
        action = policy.bayes_search(max_num_probes=1, num_search_each_probe=5,
                                     simulator=None, score='EI', interval=0, num_rand_basis=0)
    for i_action in action:
        qsub_action(filename[i_action], i_action)


main()
