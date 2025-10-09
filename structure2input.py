import pymatgen
import seekpath
import numpy
import math
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from write_openmx import write_openmx
from write_pwx import write_pwx
from write_pp import write_pp
from write_ph import write_ph
from write_wannier import write_wannier
from write_sh import write_sh
# from write_hilapw import write_hilapw
from pymatgen.core.periodic_table import get_el_sp


def structure2input(structure, dk_path, dq_grid, pseudo_kind, host, rel, uc):

    if pseudo_kind == "sg15":
        if rel:
            from sg15_rel import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
        else:
            from sg15 import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
    elif pseudo_kind == "pslibrary":
        if rel:
            from pslibrary_rel import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
        else:
            from pslibrary import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
    elif pseudo_kind == "sssp":
        from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
    elif pseudo_kind == "ssspsol":
        from ssspsol import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
    elif pseudo_kind == "sssp_us":
        from sssp_us import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
    else:
        from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict, band_dict
        print("Unsupported pseudo potential library :", pseudo_kind)
        exit(1)
    #
    # Band path and primitive lattice
    #
    frac_coord2 = numpy.array(structure.frac_coords)
    for ipos in range(len(frac_coord2)):
        for iaxis in range(3):
            coord3 = frac_coord2[ipos, iaxis] * 6.0
            if abs(round(coord3) - coord3) < 0.001:
                frac_coord2[ipos, iaxis] = float(round(coord3)) / 6.0
    #
    if uc:
        skp = seekpath.get_path((structure.lattice.matrix, frac_coord2,
                                [pymatgen.core.Element(str(spc)).number for spc in structure.species]))
        avec = skp["primitive_lattice"]
        bvec = skp["reciprocal_primitive_lattice"]
        atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
        pos = skp["primitive_positions"]
    else:
        skp = seekpath.get_path_orig_cell((structure.lattice.matrix, frac_coord2,
                                          [pymatgen.core.Element(str(spc)).number for spc in structure.species]))
        avec = structure.lattice.matrix
        bvec = structure.lattice.reciprocal_lattice.matrix
        atom = [str(get_el_sp(iat)) for iat in structure.atomic_numbers]
        pos = frac_coord2
    #
    # Lattice information
    #
    typ = sorted(set(atom))
    print("Bravais lattice : ", skp["bravais_lattice"])
    print("Space group : ", skp['spacegroup_international'])
    bz_volume = abs(- bvec[0][2]*bvec[1][1]*bvec[2][0]
                    + bvec[0][1]*bvec[1][2]*bvec[2][0] 
                    + bvec[0][2]*bvec[1][0]*bvec[2][1]
                    - bvec[0][0]*bvec[1][2]*bvec[2][1]
                    - bvec[0][1]*bvec[1][0]*bvec[2][2]
                    + bvec[0][0]*bvec[1][1]*bvec[2][2])
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
    numpw = 4.0*math.pi*math.sqrt(ecutwfc)**3/3.0/bz_volume
    print("Estimated number of PW (WFC) :", numpw)
    #
    # k and q grid
    #  n_i = c / |a_i| / dq
    #  c = (|a_1| |a_2| |a_3| / V_{uc})^{1/3}
    #
    v_uc = abs(numpy.linalg.det(avec))
    nq = numpy.zeros(3, numpy.int_)
    norm = numpy.zeros(3, numpy.float64)
    for ii in range(3):
        norm[ii] = numpy.sqrt(numpy.dot(avec[ii][:], avec[ii][:]))
    c = (norm[0]*norm[1]*norm[2] / v_uc)**(1.0/3.0)
    for ii in range(3):
        nq[ii] = round(c / norm[ii] / dq_grid)
        if nq[ii] == 0:
            nq[ii] = 1
    print("Coarse grid : ", nq[0], nq[1], nq[2])
    #
    # Get explicit kpath applicable to bands.x and plotband.x
    #
    kpath = []
    nkpath = []
    #
    dk_jump = 1.0e10
    for ipath in range(len(skp['path'])):
        dk = numpy.array(skp['point_coords'][skp['path'][ipath][1]])\
           - numpy.array(skp['point_coords'][skp['path'][ipath][0]])
        dk = numpy.dot(dk, bvec)
        dknorm = math.sqrt(numpy.dot(dk, dk))
        nkpath0 = max(2, int(dknorm / dk_path))
        nkpath.append(nkpath0)
        for ik in range(nkpath0):
            xkpath = ik / nkpath0
            kpath0 = numpy.array(skp['point_coords'][skp['path'][ipath][0]]) * (1.0 - xkpath) \
                + numpy.array(skp['point_coords'][skp['path'][ipath][1]]) * xkpath
            kpath.append(kpath0)
        #
        # jump case
        #
        if ipath < len(skp['path']) - 1 and skp['path'][ipath][1] != skp['path'][ipath+1][0]:
            kpath0 = numpy.array(skp['point_coords'][skp['path'][ipath][1]])
            kpath1 = numpy.array(skp['point_coords'][skp['path'][ipath+1][0]])
            dk = kpath1 - kpath0
            dk = numpy.dot(dk, bvec)
            dk_jump = min(dk_jump, math.sqrt(numpy.dot(dk, dk)))
            kpath.append(kpath0)
    #
    # If the jump is relatively small, shift the last point closer to the end
    #
    dk = numpy.array(kpath[1]) - numpy.array(kpath[0])
    dk = numpy.dot(dk, bvec)
    dknorm = math.sqrt(numpy.dot(dk, dk))
    if dk_jump < 10.0 * dknorm:
        xkpath = dk_jump * 0.09
        for i in range(3):
            kpath[1][i] = kpath[0][i] + (kpath[1][i] - kpath[0][i])*xkpath/dknorm
    #
    # Last point
    #
    kpath.append(numpy.array(skp['point_coords'][skp['path'][len(skp['path']) - 1][1]]))
    #
    # Band path
    #
    print("Band path")
    for ipath in range(len(skp["path"])):
        print("%5d %8s %10.5f %10.5f %10.5f %8s %10.5f %10.5f %10.5f" % (
            nkpath[ipath],
            skp['path'][ipath][0],
            skp['point_coords'][skp['path'][ipath][0]][0],
            skp['point_coords'][skp['path'][ipath][0]][1],
            skp['point_coords'][skp['path'][ipath][0]][2],
            skp['path'][ipath][1],
            skp['point_coords'][skp['path'][ipath][1]][0],
            skp['point_coords'][skp['path'][ipath][1]][1],
            skp['point_coords'][skp['path'][ipath][1]][2]))
    #
    # Number of bands
    #
    nbnd = 0
    for iat in atom:
        nbnd += valence_dict[iat]
    nbnd = nbnd + len(atom)*20
    if rel:
        nbnd *= 2
    print("Number of Bands with high-energy : ", nbnd)
    #
    # Number of bands
    #
    nbnd0 = 0
    for iat in atom:
        nbnd0 += band_dict[iat]
    print("Number of Bands : ", nbnd0)
    #
    # Shell scripts
    #
    if uc:
        structure2 = pymatgen.core.Structure(skp["primitive_lattice"],
                                             skp["primitive_types"], skp["primitive_positions"])
    else:
        structure2 = structure
    spg_analysis = SpacegroupAnalyzer(structure2)
    coarse = spg_analysis.get_ir_reciprocal_mesh(mesh=(nq[0], nq[1], nq[2]), is_shift=(0, 0, 0))
    middle = spg_analysis.get_ir_reciprocal_mesh(mesh=(nq[0]*2, nq[1]*2, nq[2]*2), is_shift=(0, 0, 0))
    dense = spg_analysis.get_ir_reciprocal_mesh(mesh=(nq[0]*4, nq[1]*4, nq[2]*4), is_shift=(0, 0, 0))
    print("Number of irreducible k : ", len(coarse), len(middle), len(dense))
    write_sh(nq[0]*nq[1]*nq[2], len(middle), len(dense),
             len(kpath), atom, atomwfc_dict, host, numpw*nbnd, rel)
    #
    # rx.in, scf.in, nscf.in, band.in , nscf_w.in, nscf_r.in
    #
    write_pwx(avec, atom, pos, ecutwfc, ecutrho, pseudo_dict, nq, nbnd0, nbnd, rel, kpath)
    #
    # ph.in, elph.in, epmat.in, phdos.in, rpa.in, scdft.in
    #
    write_ph(nq, ecutrho, nbnd)
    #
    # bands.in, pp.in, proj.in, pw2wan.in, q2r.in
    #
    write_pp()
    #
    # band.gp, pwscf.win, respack.in, disp.in
    #
    write_wannier(avec, bvec, atom, pos, skp, nbnd0, nq, atomwfc_dict, kpath)
    #
    # openmx.in : Input file for openmx
    #
    write_openmx(avec, atom, pos, skp, nq, rel, nkpath)
    #
    # HiLAPW input
    #
    # write_hilapw(skp, nq)
