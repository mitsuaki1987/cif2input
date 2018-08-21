import os
import math


def good_proc(nproc, ncore):
    if ncore == 24:
        if nproc == 5:
            nproc = 4
        elif nproc == 7:
            nproc = 6
        elif 8 < nproc < 12:
            nproc = 8
        elif 12 < nproc < 24:
            nproc = 12
    else:  # ncore == 40
        if nproc == 3:
            nproc = 2
        elif 5 < nproc < 8:
            nproc = 5
        elif nproc == 9:
            nproc = 8
        elif 10 < nproc < 20:
            nproc = 10
        elif 20 < nproc < 40:
            nproc = 20

    return nproc


def write_sh(nks, nkd, nk_path, atom, atomwfc_dict, queue):
    pw = "~/program/QE/qe-6.2.1/bin/pw.x"
    proj = "~/program/QE/qe-6.2.1/bin/projwfc.x"
    vf = "~/program/QE/qe-6.2.1/bin/fermi_velocity.x"
    bands = "~/program/QE/qe-6.2.1/bin/bands.x"
    sumpdos = "~/program/QE/qe-6.2.1/bin/sumpdos.x"
    fproj = "~/program/QE/qe-6.2.1/bin/fermi_proj.x"
    typ = set(atom)
    #
    if queue == "F4cpus":
        maxnode = 4
        ncore = 24
    elif queue == "F4cpue":
        maxnode = 4
        ncore = 40
    elif queue == "F36cpus":
        maxnode = 36
        ncore = 24
    elif queue == "F9cpue":
        maxnode = 9
        ncore = 40
    elif queue == "F36cpue":
        maxnode = 36
        ncore = 40
    else:  # queue == "F144cpus":
        maxnode = 144
        ncore = 24
    #
    # Structure optimization
    #
    nk = min(ncore*maxnode, nks)
    ntg = good_proc(int(ncore*maxnode / nk), ncore)
    nproc = nk*ntg
    node = math.ceil(nproc / ncore)
    if not os.path.isfile("rx.sh"):
        with open("rx.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in rx.in > rx_s.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("sed -n -e '/occupations/c occupations=\"tetrahedra_opt\"' -e '1,/CELL_PARAMETERS/p' rx.in > rx_t.in",
                  file=f)
            print("grep -A 3 CELL_PARAMETERS rx_s.out | tail -n 3 >> rx_t.in", file=f)
            print("awk '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS/' rx.in >> rx_t.in", file=f)
            print("grep -A %d ATOMIC_POSITIONS rx_s.out |tail -n %d >> rx_t.in" % (len(atom), len(atom)), file=f)
            print("sed -n -e '/K_POINTS/,$p' rx.in >> rx_t.in", file=f)
            print("sed -i -e '/occupations/c occupations=\"tetrahedra_opt\"' rx_t.in", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in rx_t.in > rx_t.out"
                  % (nproc, pw, nk, ntg), file=f)
    #
    # Charge optimization
    #
    if not os.path.isfile("scf.sh"):
        with open("scf.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in scf.in > scf.out"
                  % (nproc, pw, nk, ntg), file=f)
    #
    # Projected DOS
    #
    nk = min(ncore*maxnode, nkd)
    ntg = good_proc(int(ncore*maxnode / nk), ncore)
    nproc = nk * ntg
    node = math.ceil(nproc / ncore)
    #
    # Atomwfc dictionary for fermi_proj.x
    #
    pfermi = {ityp: [[] for il in range(len(atomwfc_dict[ityp][0]))] for ityp in typ}
    ii = 0
    for iat in atom:
        for il in range(len(atomwfc_dict[iat][0])):
            for im in range(atomwfc_dict[iat][0][il]):
                ii += 1
                pfermi[iat][il].append(ii)
    #
    if not os.path.isfile("proj.sh"):
        with open("proj.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in nscf.in > nscf.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("mpijob -n 1 %s -in nscf.in > vfermi.out" % vf, file=f)
            print("ef=`grep Fermi nscf.out| awk '{print $5}'`", file=f)
            print("sed -i -e '/emin/c emin = '${ef}'' -e '/emax/c emax = '${ef}'' proj.in", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in proj.in > proj.out"
                  % (nproc, proj, nk, ntg), file=f)
            #
            # Sum PDOS at each Atom and L
            #
            for ityp in typ:
                for il in range(len(atomwfc_dict[ityp][1])):
                    print("sumpdos pwscf.pdos_atm*\\(%s\\)_wfc#%d* > pdos_%s%s"
                          % (ityp, il+1, ityp, atomwfc_dict[ityp][1][il]), file=f)
            #
            # Fermi surface with atomic projection
            #
            for ityp in typ:
                for il in range(len(atomwfc_dict[ityp][1])):
                    print("sed -e '$a %d\\n" % len(pfermi[ityp][il]), end="", file=f)
                    for ii in pfermi[ityp][il]:
                        print(" %d" % ii, end="", file=f)
                    print("' proj.in > proj_f.in", file=f)
                    print("mpijob -n 1 %s -in proj_f.in" % fproj, file=f)
                    print("mv proj.frmsf %s%s.frmsf" % (ityp, atomwfc_dict[ityp][1][il]), file=f)
    #
    # Band
    #
    nk = min(ncore*maxnode, nk_path)
    ntg = good_proc(int(ncore*maxnode / nk), ncore)
    nproc = nk * ntg
    node = math.ceil(nproc / ncore)
    if not os.path.isfile("band.sh"):
        with open("band.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in band.in > band.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in bands.in > bands.out"
                  % (nproc, bands, nk_path, ntg), file=f)
