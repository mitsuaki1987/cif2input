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


def write_sh(nkcbz, nks, nkd, nk_path, atom, atomwfc_dict, host, npw_nbnd, rel):
    pw = "~/bin/pw.x"
    ph = "~/bin/ph.x"
    proj = "~/bin/projwfc.x"
    vf = "~/bin/fermi_velocity.x"
    bands = "~/bin/bands.x"
    sumpdos = "~/bin/sumpdos.x"
    fproj = "~/bin/fermi_proj.x"
    sctk = "~/bin/sctk.x"
    typ = sorted(set(atom))
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
    qsystem = ""
    extension = ""
    #
    nmem = npw_nbnd * nkcbz * 4.0**3 * 16.0 / 1.0e9 * 50.0  # Last is scale factor
    print("Estimated memory for WFC (GB) : ", nmem)
    if host == "wisteria":
        extension = "_w.sh"
        qsystem = "pj"
        mem_per_node = 32
        core_per_node = 48
        required_node = max(int(nmem / mem_per_node), 1)
        if required_node <= 72:
            queue = "short-o"
        elif required_node <= 2304:
            queue = "regular"
        else:
            print("Too large system")
            exit(-1)
        maxnode = required_node
        jobscript_queue = "#PJM -L rscgrp=" + queue
        jobscript_node = "#PJM -L node="
        jobscript_mpi = "#PJM --mpi proc="
        jobscript_omp = "#PJM --omp thread="
        jobscript_time = "#PJM -L elapse="
        jobscript_workdir = "${PJM_O_WORKDIR}"
        mpiexec = "mpiexec"
    elif host == "ohtaka":
        extension = "_o.sh"
        qsystem = "slurm"
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
    elif host == "kugui":
        extension = "_k.sh"
        qsystem = "pbs"
        mem_per_node = 256
        core_per_node = 128
        required_node = max(int(nmem / mem_per_node), 1)
        if required_node <= 4:
            queue = "F1cpu"
            maxnode = 1
        elif required_node <= 4:
            queue = "F4cpu"
            maxnode = 4
        elif required_node <= 16:
            queue = "F16cpu"
            maxnode = 16
        else:
            print("Too large system")
            exit(-1)
        jobscript_queue = "#PBS -q " + queue
        jobscript_node = "#PBS -node "
        jobscript_mpi = "#PBS -mpi "
        jobscript_omp = "#PBS -omp "
        jobscript_time = "#PBS -l walltime="
        jobscript_workdir = "${PBS_O_WORKDIR}"
        mpiexec = "mpijob"

    elif host == "imr":
        extension = "_i.sh"
        qsystem = "pbs"
        mem_per_node = 256
        core_per_node = 128
        required_node = max(int(nmem / mem_per_node), 1)
        queue = "P_030"
        maxnode = 1
        jobscript_queue = "#PBS -q " + queue
        jobscript_node = "#PBS -node "
        jobscript_mpi = "#select=1:ncpus="
        jobscript_omp = "#PBS -omp "
        jobscript_time = "#PBS -l walltime="
        jobscript_workdir = "${PBS_O_WORKDIR}"
        mpiexec = "mpirun"
    else:
        print("Unsupported host")
        exit(-1)
    print("Required node :", nmem / mem_per_node)
    #
    # Structure optimization
    #
    npool: int = min(core_per_node*maxnode, nks)
    nth = good_proc(int(core_per_node*maxnode / npool), core_per_node)
    ncore = npool*nth
    node = math.ceil(ncore / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
    with open("rx" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node/nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        print("#sed -n -e '1,/CELL_PARAMETERS/p' rx.in > temp", file=f)
        print("#grep -A 3 CELL_PARAMETERS rx.out | tail -n 3 >> temp", file=f)
        print("#awk '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS/' rx.in >> temp", file=f)
        print("#grep -A %d ATOMIC_POSITIONS rx.out |tail -n %d >> temp" % (len(atom), len(atom)), file=f)
        print("#sed -n -e '/K_POINTS/,$p' rx.in >> temp", file=f)
        print("#mv temp rx.in", file=f)
        if qsystem == "pj":
            print(mpiexec, "-of rx.out -n", npool, pw, "-npool", npool, "-in rx.in ", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in rx.in > rx.out", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
        print("find ./ -name \"*.upf\" -delete", file=f)
        print("find ./ -name \"*.UPF\" -delete", file=f)
    #
    # Charge optimization
    #
    with open("scf" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node * node / nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        if qsystem == "pj":
            print(mpiexec, "-of scf.out -n", npool, pw, "-npool", npool, "-in scf.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in scf.in > scf.out", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
    #
    # Phonon
    #
    with open("ph" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
            print(jobscript_time + "24:00:00", file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
            if queue == "short-o":
                print(jobscript_time + "8:00:00", file=f)
            else:
                print(jobscript_time + "24:00:00", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
            print(jobscript_time + "24:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        if qsystem == "pj":
            print(mpiexec, "-of nscf_p.out -n", npool, pw, "-npool", npool, "-in nscf_p.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in nscf_p.in > nscf_p.out", file=f)
        #
        # Compute number of q
        #
        print("sed -i -e \"/only_init/c only_init = .true.\" ph.in", file=f)
        print(mpiexec, "-n", npool, ph, "-npool", npool, "-in ph.in", file=f)
        print("sed -i -e \"/only_init/c only_init = .false.\" ph.in", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        #
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s matdyn${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" ph.in", file=f)
        if qsystem == "pj":
            print(" ", mpiexec, "-of ph${i}.out -n", npool, ph, "-npool", npool, "-in ph.in", file=f)
        else:
            print(" ", mpiexec, "-n", npool, ph, "-npool", npool, "-in ph.in > ph${i}.out", file=f)
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
        print("  test -s matdyn${i} || exit", file=f)
        print("done", file=f)
        print("touch PH_DONE", file=f)
    #
    # Electron-phonon
    #
    with open("elph" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd", jobscript_workdir, file=f)
        if qsystem == "pj":
            print(mpiexec, "-of nscf_pd.out -n", npool, pw, "-npool", npool, "-in nscf_pd.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in nscf_pd.in > nscf_pd.out", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s matdyn${i}.elph.${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" elph.in", file=f)
        if qsystem == "pj":
            print(" ", mpiexec, "-of elph${i}.out -n", npool, ph, "-npool", npool, "-in elph.in", file=f)
        else:
            print(" ", mpiexec, "-n", npool, ph, "-npool", npool, "-in elph.in > elph${i}.out", file=f)
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
        print("  test -s matdyn${i}.elph.${i} || exit", file=f)
        print("done", file=f)
        print("touch ELPH_DONE", file=f)
    #
    # Electron-phonon matrix
    #
    with open("epmat" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
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
        if qsystem == "pj":
            print(mpiexec, "-of nscf_pc.out -n", npool, pw, "-npool", npool, "-in nscf_pc.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in nscf_pc.in > nscf_pc.out", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s _ph0/pwscf.q_${i}/pwscf.elph${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" epmat.in", file=f)
        if qsystem == "pj":
            print(" ", mpiexec, "-of epmat${i}.out -n", npool, ph, "-npool", npool, "-in epmat.in", file=f)
        else:
            print(" ", mpiexec, "-n", npool, ph, "-npool", npool, "-in epmat.in > epmat${i}.out", file=f)
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
    npool = min(core_per_node*maxnode, nkd)
    nth = good_proc(int(core_per_node*maxnode / npool), core_per_node)
    ncore = npool * nth
    node = math.ceil(ncore / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
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
                if rel:
                    ii += 1
                    pfermi[iat][il[0]].append(ii)
    #
    with open("proj" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        if qsystem == "pj":
            print(mpiexec, "-of nscf.out -n", npool, pw, "-npool", npool, "-in nscf.in", file=f)
            print(mpiexec, "-of vfermi.out -n 1", vf, "-in nscf.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in nscf.in > nscf.out", file=f)
            print(mpiexec, "-n 1", vf, "-in nscf.in > vfermi.out", file=f)
        print("ef=`grep Fermi nscf.out| awk '{print $5}'`", file=f)
        print("sed -i -e '/emin/c emin = '${ef}'' -e '/emax/c emax = '${ef}'' proj.in", file=f)
        if qsystem == "pj":
            print(mpiexec, "-of proj.out -n", npool, proj, "-npool", npool, "-in proj.in", file=f)
        else:
            print(mpiexec, "-n", npool, proj, "-npool", npool, "-in proj.in > proj.out", file=f)
        #
        # Sum PDOS at each Atom and L
        #
        for ityp in typ:
            nwfc = 1
            for il in atomwfc_dict[ityp]:
                if rel and int(il[1]) > 1:
                    print("%s pwscf.pdos_atm*\\(%s\\)_wfc#%d* pwscf.pdos_atm*\\(%s\\)_wfc#%d* > pdos_%s%s"
                          % (sumpdos, ityp, nwfc, ityp, nwfc+1, ityp, il[0]), file=f)
                    nwfc += 2
                else:
                    print("%s pwscf.pdos_atm*\\(%s\\)_wfc#%d* > pdos_%s%s"
                          % (sumpdos, ityp, nwfc, ityp, il[0]), file=f)
                    nwfc += 1
        #
        # Fermi surface with atomic projection
        #
        for ityp in typ:
            for il in atomwfc_dict[ityp]:
                print("sed -e '$a %d\\n" % len(pfermi[ityp][il[0]]), end="", file=f)
                for ii in pfermi[ityp][il[0]]:
                    print(" %d" % ii, end="", file=f)
                print("' proj.in > proj_f.in", file=f)
                print(mpiexec, "-n 1", fproj, "-in proj_f.in", file=f)
                print("mv pwscf_proj.frmsf " + ityp + il[0] + ".frmsf", file=f)
        print("find ./ -name \"pwscf.wfc*\" -delete", file=f)
        print("find ./ -name \"wfc*.dat\" -delete", file=f)
    #
    # Band
    #
    npool = min(core_per_node*maxnode, nk_path)
    nth = good_proc(int(core_per_node*maxnode / npool), core_per_node)
    ncore = npool * nth
    node = math.ceil(ncore / core_per_node)
    if host == "ohtaka" and node > 36:
        node = maxnode
    with open("band" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
        print(jobscript_time + "8:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        if qsystem == "pj":
            print(mpiexec, "-of band.out -n", npool, pw, "-npool", nk_path, "-in band.in", file=f)
            print(mpiexec, "-of bandsx.out -n", npool, bands, "-npool", nk_path, "-in bands.in", file=f)
            print(mpiexec, "-of pband.out -n", npool, proj, "-npool", npool, "-in proj.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", nk_path, "-in band.in > band.out", file=f)
            print(mpiexec, "-n", npool, bands, "-npool", nk_path, "-in bands.in > bandsx.out", file=f)
            print(mpiexec, "-n", npool, proj, "-npool", npool, "-in proj.in > pband.out", file=f)
        #
        # Projected band
        #
        print("mv bands.projwfc_up bands.out.proj", file=f)
        for ityp in typ:
            for il in atomwfc_dict[ityp]:
                print("sed -e '2i", end="", file=f)
                for ii in pfermi[ityp][il[0]]:
                    print(" %d" % ii, end="", file=f)
                print("' -e '3c " + ityp + il[0] + ".xmgr' plotband.in > plotpband.in", file=f)
                if qsystem == "pj":
                    print(mpiexec, "-n 1 -stdin plotpband.in ~/bin/plotband.x", file=f)
                else:
                    print(mpiexec, "-n 1 ~/bin/plotband.x < plotpband.in", file=f)
    #
    # Coulomb matrix
    #
    node = maxnode
    nth = core_per_node / 4
    npool = maxnode * 4
    if host == "ohtaka" and node > 36:
        node = maxnode
    with open("kel" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            proc_per_node = int(core_per_node / nth)
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=" + str(proc_per_node) + ":ompthreads=" + str(nth), file=f)
            print(jobscript_time + "24:00:00", file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
            if queue == "short-o":
                print(jobscript_time + "8:00:00", file=f)
            else:
                print(jobscript_time + "24:00:00", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_omp, nth, file=f)
            print(jobscript_time + "24:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        if qsystem == "pj":
            print(mpiexec, "-of twin.out -n", npool, pw, "-npool", npool, "-in twin.in", file=f)
        else:
            print(mpiexec, "-n", npool, pw, "-npool", npool, "-in twin.in > twin.out", file=f)
        print("sed -i -e \'/calculation/c calculation = \"kel\"\' sctk.in", file=f)
        print("nq=`awk \'NR==2{print $1}\' matdyn0`", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s pwscf.vel${i} && continue", file=f)
        print("  sed -i -e \"/start_q/c start_q=$i\" -e \"/last_q/c last_q=$i\" sctk.in", file=f)
        if qsystem == "pj":
            print(" ", mpiexec, "-of kel${i}.out -n", npool, sctk, "-npool", npool, "-in sctk.in", file=f)
        else:
            print(" ", mpiexec, "-n", npool, sctk, "-npool", npool, "-in sctk.in > kel${i}.out", file=f)
        print("done", file=f)
        print("for i in `seq 1 ${nq}`; do", file=f)
        print("  test -s pwscf.vel${i} || exit", file=f)
        print("done", file=f)
        print("touch KEL_DONE", file=f)
    #
    # 0K calculation, lambda_mu_k, delta_f
    #
    node = maxnode
    with open("scdft0" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=1:ompthreads=" + str(core_per_node), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_mpi, node, file=f)
            print(jobscript_omp, core_per_node, file=f)
        print(jobscript_time + "3:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print("sed -i -e \'/calculation/c calculation = \"lambda_mu_k\"\' sctk.in", file=f)
        if qsystem == "pj":
            print(mpiexec, "-of lambda_mu_k.out -n", node, sctk, "-npool", node, "-in sctk.in", file=f)
        else:
            print(mpiexec, "-n", node, sctk, "-npool", node, "-in sctk.in > lambda_mu_k.out", file=f)
        print("sed -i -e \'/calculation/c calculation = \"scdft\"\' sctk.in", file=f)
        if qsystem == "pj":
            print(mpiexec, "-of scdft0.out -n", node, sctk, "-npool", node, "-in sctk.in", file=f)
        else:
            print(mpiexec, "-n", node, sctk, "-npool", node, "-in sctk.in > scdft0.out", file=f)
        print("sed -i -e \'/calculation/c calculation = \"deltaf\"\' sctk.in", file=f)
        if qsystem == "pj":
            print(mpiexec, "-of deltaf.out -n", node, sctk, "-npool", node, "-in sctk.in", file=f)
        else:
            print(mpiexec, "-n", node, sctk, "-npool", node, "-in sctk.in > deltaf.out", file=f)
    #
    # Tc calculation
    #
    node = maxnode
    with open("scdft" + extension, 'w') as f:
        print("#!/bin/sh", file=f)
        print(jobscript_queue, file=f)
        if qsystem == "pbs":
            print("#PBS -l select=" + str(node) + ":ncpus=" + str(core_per_node)
                  + ":mpiprocs=1:ompthreads=" + str(core_per_node), file=f)
        elif qsystem == "pj":
            print(jobscript_node + str(node), file=f)
            print(jobscript_mpi + str(int(core_per_node*node/nth)), file=f)
            print(jobscript_omp + str(nth), file=f)
            print("#PJM -g ga20", file=f)
            print("#PJM -j", file=f)
            print("#PJM -x PJM_FEFS_CACHE_MODE=3", file=f)
        else:
            print(jobscript_node, node, file=f)
            print(jobscript_mpi, node, file=f)
            print(jobscript_omp, core_per_node, file=f)
        print(jobscript_time + "3:00:00", file=f)
        print("source ~/.bashrc", file=f)
        print("cd ", jobscript_workdir, file=f)
        print("sed -i -e \'/calculation/c calculation = \"scdft_tc\"\' sctk.in", file=f)
        if qsystem == "pj":
            print(mpiexec, "-of tc.out -n", node, sctk, "-npool", node, "-in sctk.in", file=f)
        else:
            print(mpiexec, "-n", node, sctk, "-npool", node, "-in sctk.in > tc.out", file=f)
