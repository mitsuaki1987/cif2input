import os


def write_sh(nks, nkd, nk_path, atom, prefix, atomwfc_dict):
    pw = "~/program/QE/qe-6.2.1/bin/pw.x"
    proj = "~/program/QE/qe-6.2.1/bin/projwfc.x"
    vf = "~/program/QE/qe-6.2.1/bin/fermi_velocity.x"
    bands = "~/program/QE/qe-6.2.1/bin/bands.x"
    sumpdos = "~/program/QE/qe-6.2.1/bin/sumpdos.x"
    fproj = "~/program/QE/qe-6.2.1/bin/fermi_proj.x"
    typ = set(atom)
    #
    # Structure optimization
    #
    nk = min(24*4, nks)
    ntg = int(24*4 / nk)  # 2, 3, 4, 6, 8, 12, 24
    if ntg == 5:
        ntg = 4
    elif ntg == 7:
        ntg = 6
    elif 8 < ntg < 12:
        ntg = 8
    elif 12 < ntg < 24:
        ntg = 12
    nproc = nk*ntg
    node = int(nproc / 24)
    if node <= 4:
        queue = "F4cpi"
    elif node <= 36:
        queue = "F36cpu"
    else:
        queue = "F144cpu"
    if not os.path.isfile("rx.sh"):
        with open("rx.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue, file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in rx.in > rx_s.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("sed -n -e '/occupations/c occupations=\"tetrahedra_opt\"' -e '1,/CELL_PARAMETERS/p' rx.in > rx_t.in",
                  file=f)
            print("awk '/Begin final coordinates/,/End final coordinate/' rx_s.out | sed -n -e '6,8p' >> rx_t.in",
                  file=f)
            print("awk '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS/' rx.in >> rx_t.in", file=f)
            print("awk '/Begin final coordinates/,/End final coordinate/' "
                  "rx_s.out | sed -e '$d'| sed -n -e '11,$p' >> rx_t.in", file=f)
            print("sed -n -e '/K_POINTS/,$p' rx.in >> rx_t.in", file=f)
            print("sed -i -e '/occupations/c occupations=\"tetrahedra_opt\"' rx_t.in", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in rx_t.in > rx_t.out"
                  % (nproc, pw, nks, ntg), file=f)
    #
    # Charge optimization
    #
    if not os.path.isfile("scf.sh"):
        with open("scf.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue, file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in scf.in > scf.out"
                  % (nproc, pw, nks, ntg), file=f)
    #
    # Projected DOS
    #
    nk = min(24*4, nkd)
    ntg = int(24*4 / nk)
    if ntg == 5:
        ntg = 4
    elif ntg == 7:
        ntg = 6
    elif 8 < ntg < 12:
        ntg = 8
    elif 12 < ntg < 24:
        ntg = 12
    nproc = nk * ntg
    node = int(nproc / 24)
    if node <= 4:
        queue = "F4cpi"
    elif node <= 36:
        queue = "F36cpu"
    else:
        queue = "F144cpu"
    #
    # Atomwfc dictionary for fermi_proj.x
    #
    pfermi = {ityp: [[] for il in range(len(atomwfc_dict[ityp]))] for ityp in typ}
    ii = 0
    for iat in atom:
        for il in range(len(atomwfc_dict[iat])):
            for im in range(atomwfc_dict[iat][il]):
                ii += 1
                pfermi[iat][il].append(ii)
    #
    if not os.path.isfile("proj.sh"):
        with open("proj.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue, file=f)
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
                for il in range(len(atomwfc_dict[ityp])):
                    print("%s %s.pdos_atm*\\(%s\\)_wfc#%d* > %s.pdos_%s%d"
                          % (sumpdos, prefix, ityp, il+1, prefix, ityp, il+1), file=f)
            #
            # Fermi surface with atomic projection
            #
            for ityp in typ:
                for il in range(len(atomwfc_dict[ityp])):
                    print("sed -e '$a %d\\n" % len(pfermi[ityp][il]), end="", file=f)
                    for ii in pfermi[ityp][il]:
                        print(" %d" % ii, end="", file=f)
                    print("' proj.in > proj_f.in", file=f)
                    print("mpijob -n 1 %s -in proj_f.in" % fproj, file=f)
                    print("mv proj.frmsf %s%d.frmsf" % (ityp, il+1), file=f)
    #
    # Band
    #
    nk = min(24*4, nk_path)
    ntg = int(24*4 / nk)
    if ntg == 5:
        ntg = 4
    elif ntg == 7:
        ntg = 6
    elif 8 < ntg < 12:
        ntg = 8
    elif 12 < ntg < 24:
        ntg = 12
    nproc = nk * ntg
    node = int(nproc / 24)
    if node <= 4:
        queue = "F4cpi"
    elif node <= 36:
        queue = "F36cpu"
    else:
        queue = "F144cpu"
    if not os.path.isfile("band.sh"):
        with open("band.sh", 'w') as f:
            print("#!/bin/sh", file=f)
            print("#QSUB -queue", queue, file=f)
            print("#QSUB -node", node, file=f)
            print("#PBS -l walltime=8:00:00", file=f)
            print("source ~/.bashrc", file=f)
            print("cd $PBS_O_WORKDIR", file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in band.in > band.out"
                  % (nproc, pw, nk, ntg), file=f)
            print("mpijob -n %d %s -nk %d -ntg %d -in bands.in > bands.out"
                  % (nproc, bands, nk_path, ntg), file=f)
