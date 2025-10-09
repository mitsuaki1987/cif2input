#!/usr/bin/python3
import os
import pymatgen
import seekpath
from pymatgen.core.periodic_table import get_el_sp
from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict
from xml.etree import ElementTree
import physbo
import subprocess
import numpy
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def load_descriptor():
    with open("desc.dat", "r") as f:
        lines = f.readlines()
        filename = []
        descriptor = []
        for line in lines:
            filename.append(str(line.split()[0]))
            descriptor.append(numpy.array(line.split()[1:], dtype=numpy.float64))

    return numpy.array(descriptor), filename


def load_result():

    with open("dos.dat", "r") as f:
        lines = f.readlines()
        action = []
        result = []
        for line in lines:
            action.append(int(line.split()[0]))
            result.append(float(line.split()[1]))

    return numpy.array(action), numpy.array(result)


class Simulator:
    def __init__(self):
        _, self.filename = load_descriptor()

    def __call__(self, action):

        structure = pymatgen.core.Structure.from_file(self.filename[action[0]])
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
        typ = sorted(set(atom))
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
        # SCF file
        #
        with open("scf" + str(action[0]) + ".in", 'w') as f:
            print("&CONTROL", file=f)
            print(" calculation = \'scf\'", file=f)
            print(" prefix = \'%s\'" % str(action[0]), file=f)
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
        #
        # Run DFT
        #
        subprocess.call("mpirun -hostfile $PBS_NODEFILE -np %d ~/bin/pw.x -nk %d -in scf%d.in > scf%d.out"
                        % (n_proc, n_proc, action[0], action[0]), shell=True)
        #
        # Extract DOS
        #
        xmlfile = os.path.join(str(action[0]) + ".save/", 'data-file-schema.xml')
        tree = ElementTree.parse(xmlfile)
        root = tree.getroot()
        child = root.find('output').find('band_structure')
        efermi = float(child.find('fermi_energy').text) * 13.60569228 * 2.0
        #
        # DOS file
        #
        with open("dos" + str(action[0]) + ".in", 'w') as f:
            print("&DOS", file=f)
            print("    prefix = \'%s\'" % str(action[0]), file=f)
            print("      emin = %f" % efermi, file=f)
            print("      emax = %f" % efermi, file=f)
            print("    deltae = 0.1", file=f)
            print("    bz_sum = \"tetrahedra_opt\"", file=f)
            print("/", file=f)
        #
        # Run DOS
        #
        subprocess.call("mpirun -np 1 ~/bin/dos.x -in dos%d.in > dos%d.out" % (action[0], action[0]), shell=True)
        #
        with open(str(action[0]) + ".dos", "r") as f:
            f.readline()
            line = f.readline()
            dos = float(line.split()[1]) / float(nat)
        #
        with open("dos.dat", "a") as f:
            print(action[0], dos, self.filename[action[0]], file=f)
        #
        subprocess.call("rm -rf %d.save" % action, shell=True)
        #
        return dos


def main():
    descriptor, filename = load_descriptor()
    descriptor = physbo.misc.centering(descriptor)
    policy = physbo.search.discrete.policy(test_X=descriptor)
    policy.set_seed(1)
    #
    # Read previous result
    #
    action, result = load_result()
    if len(action) == 0:
        random_search = policy.random_search(max_num_probes=5, simulator=Simulator())
    else:
        policy.write(action, result)
    bayes_search = policy.bayes_search(max_num_probes=30, simulator=Simulator(), score='EI',
                                       interval=1, num_rand_basis=5000)

    print('f(x)=')
    print(bayes_search.fx[0:bayes_search.total_num_search])
    best_fx, best_action = bayes_search.export_all_sequence_best_fx()
    print('current best')
    print(best_fx)
    print('current best action=')
    print(best_action)
    print('history of chosen actions=')
    print(bayes_search.chosed_actions[0:bayes_search.total_num_search])


main()
