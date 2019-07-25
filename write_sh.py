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


def write_sh(nkcbz, nkc, nks, nkd, nk_path, atom, atomwfc_dict, queue):
    pw = "~/bin/pw.x"
    ph = "~/bin/ph.x"
    proj = "~/bin/projwfc.x"
    vf = "~/bin/fermi_velocity.x"
    bands = "~/bin/bands.x"
    sumpdos = "~/bin/sumpdos.x"
    fproj = "~/bin/fermi_proj.x"
    sctk = "~/bin/sctk.x"
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
    # Phonon
    #
    if not os.path.isfile("ph.sh"):
        with open("ph.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in nscf_p.in > nscf_p.out"
                  % (nproc, pw, nk, ntg), file=f)
            #
            # Compute number of q
            #
            print("sed -i -e \"/only_init/c only_init = .true.\" ph.in", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in ph.in >> ph.out"
                  % (nproc, ph, nk, ntg), file=f)
            print("sed -i -e \"/only_init/c only_init = .false.\" ph.in", file=f)
            print("nq=`awk \'NR==2{print $1}\' matdyn0`" % nkcbz, file=f)
            #
            print("for i in `seq 1 ${nq}`; do" % nkcbz, file=f)
            print("  test -s matdyn${i} && continue", file=f)
            print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" ph.in", file=f)
            print("  mpijob -n %d %s -nk %d -ntg %d -in ph.in >> ph.out"
                  % (nproc, ph, nk, ntg), file=f)
            print("  find ./ -name \"*.wfc*\" -delete", file=f)
            print("done", file=f)
            print("touch PH_DONE", file=f)
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
                    print("%s pwscf.pdos_atm*\\(%s\\)_wfc#%d* > pdos_%s%s"
                          % (sumpdos, ityp, il+1, ityp, atomwfc_dict[ityp][1][il]), file=f)
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
    # Electron-phonon
    #
    if not os.path.isfile("elph.sh"):
        with open("elph.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in nscf_pd.in > nscf_pd.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("nq=`awk \'NR==2{print $1}\' matdyn0`" % nkcbz, file=f)
            print("for i in `seq 1 ${nq}`; do" % nkcbz, file=f)
            print("  test -s matdyn${i}.elph.${i} && continue", file=f)
            print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" elph.in", file=f)
            print("  mpijob -n %d %s -nk %d -ntg %d -in elph.in >> elph.out"
                  % (nproc, ph, nk, ntg), file=f)
            print("  find ./ -name \"*.wfc*\" -delete", file=f)
            print("done", file=f)
            print("touch ELPH_DONE", file=f)
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
    #
    # Electron-phonon matrix
    #
    nk = min(ncore * maxnode, nkc)
    ntg = int(good_proc(int(ncore * maxnode / nk), ncore) / 2)
    nproc = nk * ntg
    node = math.ceil(nproc / ncore)
    if not os.path.isfile("epmat.sh"):
        with open("epmat.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("bmax=`grep \"Highest band which contains FS\" vfermi.out elph.out| awk 'NR==1{print $NF}'`",
                  file=f)
            print("bmin=`grep \"Lowest band which contains FS\" vfermi.out elph.out| awk 'NR==1{print $NF}'`",
                  file=f)
            print("sed -i -e \"/elph_nbnd_min/c elph_nbnd_min=$bmin\" "
                  "-e \"/elph_nbnd_max/c elph_nbnd_max=$bmax\" epmat.in", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in nscf_pc.in > nscf_pc.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("nq=`awk \'NR==2{print $1}\' matdyn0`" % nkcbz, file=f)
            print("for i in `seq 1 ${nq}`; do" % nkcbz, file=f)
            print("  test -s elph${i}.dat && continue", file=f)
            print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" epmat.in", file=f)
            print("  mpijob -n %d %s -nk %d -ntg %d -in epmat.in >> epmat.out"
                  % (nproc, ph, nk, ntg), file=f)
            print("  find ./ -name \"*.wfc*\" -delete", file=f)
            print("done", file=f)
            print("touch EPMAT_DONE", file=f)
    #
    # Coulomb matrix
    #
    nk = min(ncore * maxnode, nkcbz)
    ntg = good_proc(int(ncore*maxnode / nk), ncore)
    nproc = nk * ntg
    node = math.ceil(nproc / ncore)
    if not os.path.isfile("kel.sh"):
        with open("kel.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in twin.in > twin.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("export OMP_NUM_THREADS=%d" % ntg, file=f)
            print("sed -i -e \'/calculation/c calculation = \"kel\"\' sctk.in", file=f)
            print("nq=`awk \'NR==2{print $1}\' matdyn0`" % nkcbz, file=f)
            print("for i in `seq 1 ${nq}`; do" % nkcbz, file=f)
            print("  test -s vel${i}.dat && continue", file=f)
            print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" sctk.in", file=f)
            print("  mpijob -n %d %s -nk %d -in sctk.in >> kel.out"
                  % (nproc, sctk, nk), file=f)
            print("done", file=f)
            print("touch KEL_DONE", file=f)
    #
    # Coulomb matrix
    #
    node = maxnode
    if not os.path.isfile("scdft0.sh"):
        with open("scdft0.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#QSUB -mpi", node, file=f)
            print("#QSUB -omp", ncore, file=f)
            print("#PBS -l walltime=3:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("sed -i -e \'/calculation/c calculation = \"lambda_mu_k\"\' sctk.in", file=f)
            print("mpijob -n %d %s -nk %d -in sctk.in > lambda_mu_k.out"
                  % (node, sctk, node), file=f)
            print("sed -i -e \'/calculation/c calculation = \"scdft\"\' sctk.in", file=f)
            print("mpijob -n %d %s -nk %d -in sctk.in > scdft0.out"
                  % (node, sctk, node), file=f)
    #
    # Coulomb matrix
    #
    node = maxnode
    if not os.path.isfile("scdft.sh"):
        with open("scdft.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue[0:len(queue) - 1], file=f)
            print("#QSUB -node", node, file=f)
            print("#QSUB -mpi", node, file=f)
            print("#QSUB -omp", ncore, file=f)
            print("#PBS -l walltime=3:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("sed -i -e \'/calculation/c calculation = \"scdft_tc\"\' sctk.in", file=f)
            print("mpijob -n %d %s -nk %d -in sctk.in > tc.out"
                  % (node, sctk, node), file=f)
            print("sed -i -e \'/calculation/c calculation = \"deltaf\"\' sctk.in", file=f)
            print("mpijob -n %d %s -nk %d -in sctk.in > deltaf.out"
                  % (node, sctk, node), file=f)
