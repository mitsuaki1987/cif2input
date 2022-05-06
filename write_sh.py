import math


def good_proc(nproc, ncore):
    if nproc > ncore:
        nproc = ncore
    else:
        for iproc in range(nproc):
            if ncore % nproc == 0:
                break
            else:
                nproc -= 1

    return nproc


def write_sh(nkcbz, nkc, nks, nkd, nk_path, atom, atomwfc_dict, host, npw_nbnd):
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
    core_per_node = 0
    maxnode = 0
    mem_per_node = 0
    jobscript_queue = ""
    jobscript_node = ""
    jobscript_mpi = ""
    jobscript_omp = ""
    queue = ""
    jobscript_time = ""
    jobscript_workdir = ""
    mpiexec = ""
    #
    nmem = npw_nbnd * nkcbz * 4.0**3 * 16.0 / 1.0e9 * 50.0  # Last is scale factor
    print("Estimated memory for WFC (GB) : ", nmem)
    if host == "ohtaka":
        mem_per_node = 256
        core_per_node = 128
        required_node = max(int(nmem / mem_per_node), 1)
        if required_node == 1:
            queue = "F1cpu"
            maxnode = 1
        elif required_node <= 4:
            queue = "F4cpu"
            maxnode = 4
        elif required_node <= 16:
            queue = "F16cpu"
            maxnode = 16
        elif required_node <= 36:
            queue = "F36cpu"
            maxnode = 36
        elif required_node <= 72:
            queue = "F72cpu"
            maxnode = 72
        elif required_node <= 144:
            queue = "F144cpu"
            maxnode = 144
        else:
            print("Too large system")
            exit(-1)
        jobscript_queue = "#SBATCH -p " + queue
        jobscript_node = "#SBATCH -N "
        jobscript_mpi = "#SBATCH -n "
        jobscript_omp = "#SBATCH -c "
        jobscript_time = "#SBATCH -t "
        jobscript_workdir = "${SLURM_SUBMIT_DIR}"
        mpiexec = "srun"
    elif host == "enaga":
        mem_per_node = 192
        core_per_node = 40
        required_node = max(int(nmem / mem_per_node), 1)
        if required_node <= 4:
            queue = "F4cpu"
            maxnode = 4
        elif required_node <= 9:
            queue = "F9cpu"
            maxnode = 9
        elif required_node <= 36:
            queue = "F36cpu"
            maxnode = 36
        else:
            print("Too large system")
            exit(-1)
        jobscript_queue = "#QSUB -queue " + queue
        jobscript_node = "#QSUB -node "
        jobscript_mpi = "#QSUB -mpi "
        jobscript_omp = "#QSUB -omp "
        jobscript_time = "#PBS -l walltime="
        jobscript_workdir = "${PBS_O_WORKDIR}"
        mpiexec = "mpijob"
    else:
        print("Unsupported host")
        exit(-1)
    print("Required node :", nmem / mem_per_node)
    #
    # Structure optimization
    #
    nk = min(core_per_node*maxnode, nks)
    ntg = good_proc(int(core_per_node*maxnode / nk), core_per_node)
    nproc = nk*ntg
    node = math.ceil(nproc / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
    with open("rx.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        print("#sed -n -e '1,/CELL_PARAMETERS/p' rx.in > temp", file=f)
        print("#grep -A 3 CELL_PARAMETERS rx.out | tail -n 3 >> temp", file=f)
        print("#awk '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS/' rx.in >> temp", file=f)
        print("#grep -A %d ATOMIC_POSITIONS rx.out |tail -n %d >> temp" % (len(atom), len(atom)), file=f)
        print("#sed -n -e '/K_POINTS/,$p' rx.in >> temp", file=f)
        print("#mv temp rx.in", file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in rx.in > rx.out", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
        print("find ./ -name \"*.upf\" -delete", file=f)
        print("find ./ -name \"*.UPF\" -delete", file=f)
    #
    # Charge optimization
    #
    with open("scf.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in scf.in > scf.out", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
    #
    # Phonon
    #
    with open("ph.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "24:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in nscf_p.in > nscf_p.out", file=f)
        #
        # Compute number of q
        #
        print("sed -i -e \"/only_init/c only_init = .true.\" ph.in", file=f)
        print(mpiexec, "-n", nproc, ph, "-nk", nk, "-pd T -ntg", ntg, "-in ph.in", file=f)
        print("sed -i -e \"/only_init/c only_init = .false.\" ph.in", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        #
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s matdyn${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" ph.in", file=f)
        print(" ", mpiexec, "-n", nproc, ph, "-nk", nk, "-pd T -ntg", ntg, "-in ph.in >> ph.out", file=f)
        print("  find ./_ph0/ -name \"pwscf.wfc*\" -delete", file=f)
        print("  find ./_ph0/ -name \"wfc*.dat\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.dwf*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.bar*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.mix*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.recover*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.prd*\" - delete", file=f)
        print("done", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s matdyn${i} || exit", file=f)
        print("done", file=f)
        print("touch PH_DONE", file=f)
    #
    # Electron-phonon
    #
    with open("elph.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in nscf_pd.in > nscf_pd.out", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s matdyn${i}.elph.${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" elph.in", file=f)
        print(" ", mpiexec, "-n", nproc, ph, "-nk", nk, "-pd T -ntg", ntg, "-in elph.in >> elph.out", file=f)
        print("  find ./_ph0/ -name \"pwscf.wfc*\" -delete", file=f)
        print("  find ./_ph0/ -name \"wfc*.dat\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.dwf*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.bar*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.mix*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.recover*\" - delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.prd*\" - delete", file=f)
        print("done", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s matdyn${i}.elph.${i} || exit", file=f)
        print("done", file=f)
        print("touch ELPH_DONE", file=f)
    #
    # Electron-phonon matrix
    #
    with open("epmat.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        print("bmax=`grep \"Highest band which contains FS\" vfermi.out elph.out| awk 'NR==1{print $NF}'`",
              file=f)
        print("bmin=`grep \"Lowest band which contains FS\" vfermi.out elph.out| awk 'NR==1{print $NF}'`",
              file=f)
        print("sed -i -e \"/elph_nbnd_min/c elph_nbnd_min=$bmin\" "
              "-e \"/elph_nbnd_max/c elph_nbnd_max=$bmax\" epmat.in", file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in nscf_pc.in > nscf_pc.out", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s _ph0/pwscf.q_${i}/pwscf.elph${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" epmat.in", file=f)
        print(" ", mpiexec, "-n", nproc, ph, "-nk", nk, "-pd T -ntg", ntg, "-in epmat.in >> epmat.out", file=f)
        print("  find ./_ph0/ -name \"pwscf.wfc*\" -delete", file=f)
        print("  find ./_ph0/ -name \"wfc*.dat\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.dwf*\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.bar*\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.mix*\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.recover*\" -delete", file=f)
        print("  find ./_ph0/ -name \"pwscf.prd*\" -delete", file=f)
        print("done", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s _ph0/pwscf.q_${i}/pwscf.elph${i} || exit", file=f)
        print("done", file=f)
        print("touch EPMAT_DONE", file=f)
    #
    # Projected DOS
    #
    nk = min(core_per_node*maxnode, nkd)
    ntg = good_proc(int(core_per_node*maxnode / nk), core_per_node)
    nproc = nk * ntg
    node = math.ceil(nproc / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
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
    with open("proj.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in nscf.in > nscf.out", file=f)
        print(mpiexec, "-n 1", vf, "-in nscf.in > vfermi.out", file=f)
        print("ef=`grep Fermi nscf.out| awk '{print $5}'`", file=f)
        print("sed -i -e '/emin/c emin = '${ef}'' -e '/emax/c emax = '${ef}'' proj.in", file=f)
        print(mpiexec, "-n", nproc, proj, "-nk", nk, "-pd T -ntg", ntg, "-in proj.in > proj.out", file=f)
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
                print(mpiexec, "-n 1", fproj, "-in proj_f.in", file=f)
                print("mv pwscf_proj.frmsf " + ityp + atomwfc_dict[ityp][1][il] + ".frmsf", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
    #
    # Band
    #
    nk = min(core_per_node*maxnode, nk_path)
    ntg = good_proc(int(core_per_node*maxnode / nk), core_per_node)
    nproc = nk * ntg
    node = math.ceil(nproc / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
    with open("band.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk_path, "-pd T -ntg", ntg, "-in band.in > band.out", file=f)
        print(mpiexec, "-n", nproc, bands, "-nk", nk_path, "-pd T -ntg", ntg, "-in bands.in", file=f)
        print(mpiexec, "-n", nproc, proj, "-nk", nk, "-pd T -ntg", ntg, "-in proj.in > pband.out", file=f)
        #
        # Projected band
        #
        print("mv bands.projwfc_up bands.out.proj", file=f)
        for ityp in typ:
            for il in range(len(atomwfc_dict[ityp][1])):
                print("sed -e '2c", end="", file=f)
                for ii in pfermi[ityp][il]:
                    print(" %d" % ii, end="", file=f)
                print("' -e '4c " + ityp + atomwfc_dict[ityp][1][il] + ".xmgr' plotband.in > plotpband.in", file=f)
                print(mpiexec, "-n 1 ~/bin/plotband.x < plotpband.in", file=f)
    #
    # Coulomb matrix
    #
    nk = min(core_per_node * maxnode, nkcbz)
    ntg = good_proc(int(core_per_node*maxnode / nk), core_per_node)
    nproc = nk * ntg
    node = math.ceil(nproc / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
    with open("kel.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_omp, 1, file=f)
        print(jobscript_time + "24:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print(mpiexec, "-n", nproc, pw, "-nk", nk, "-pd T -ntg", ntg, "-in twin.in > twin.out", file=f)
        print("export OMP_NUM_THREADS=" + str(ntg), file=f)
        print("sed -i -e \'/calculation/c calculation = \"kel\"\' sctk.in", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s vel${i}.dat && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" sctk.in", file=f)
        print(" ", mpiexec, "-n", nproc, sctk, "-nk", nk, "-in sctk.in >> kel.out", file=f)
        print("done", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s pwscf.vel${i} || exit", file=f)
        print("done", file=f)
        print("touch KEL_DONE", file=f)
    #
    # 0K calculation, lambda_mu_k, delta_f
    #
    node = maxnode
    with open("scdft0.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_mpi, node, file=f)
        print(jobscript_omp, core_per_node, file=f)
        print(jobscript_time + "3:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print("sed -i -e \'/calculation/c calculation = \"lambda_mu_k\"\' sctk.in", file=f)
        print(mpiexec, "-n", node, sctk, "-nk", node, "-in sctk.in > lambda_mu_k.out", file=f)
        print("sed -i -e \'/calculation/c calculation = \"scdft\"\' sctk.in", file=f)
        print(mpiexec, "-n", node, sctk, "-nk", node, "-in sctk.in > scdft0.out", file=f)
        print("sed -i -e \'/calculation/c calculation = \"deltaf\"\' sctk.in", file=f)
        print(mpiexec, "-n", node, sctk, "-nk", node, "-in sctk.in > deltaf.out", file=f)
    #
    # Tc calculation
    #
    node = maxnode
    with open("scdft.sh", 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        print(jobscript_node, node, file=f)
        print(jobscript_mpi, node, file=f)
        print(jobscript_omp, core_per_node, file=f)
        print(jobscript_time + "3:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print("sed -i -e \'/calculation/c calculation = \"scdft_tc\"\' sctk.in", file=f)
        print(mpiexec, "-n", node, sctk, "-nk", node, "-in sctk.in > tc.out", file=f)
