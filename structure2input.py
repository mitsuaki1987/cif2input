import pymatgen
import seekpath
import numpy
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from write_openmx import write_openmx
from write_pwx import write_pwx
from write_pp import write_pp
from write_ph import write_ph
from write_wannier import write_wannier
from write_sh import write_sh
from pymatgen.core.periodic_table import get_el_sp


def structure2input(structure, prefix, dk_path, dq_grid, pseudo_kind, pseudo_dir):

    if pseudo_kind == "sg15":
        from sg15 import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict
    elif pseudo_kind == "pslibrary":
        from pslibrary import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict
    else:
        from sssp import pseudo_dict, ecutwfc_dict, ecutrho_dict, valence_dict, atomwfc_dict
    #
    # Band path and primitive lattice
    #
    skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                                        [pymatgen.Element(str(spc)).number for spc in structure.species]),
                                       reference_distance=dk_path)
    #
    # Lattice information
    #
    bvec = skp["reciprocal_primitive_lattice"]
    atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    typ = set(atom)
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
    # k and q grid
    #
    nq = numpy.zeros(3, numpy.int_)
    for ii in range(3):
        norm = numpy.sqrt(numpy.dot(bvec[ii][:], bvec[ii][:]))
        nq[ii] = round(norm / dq_grid)
        print(norm)
    print("Coarse grid : ", nq[0], nq[1], nq[2])
    #
    # Band path
    #
    print("Band path")
    for ipath in range(len(skp["path"])):
        start = skp["explicit_segments"][ipath][0]
        final = skp["explicit_segments"][ipath][1] - 1
        print("%5d %8s %10.5f %10.5f %10.5f %8s %10.5f %10.5f %10.5f" % (
            final - start + 1,
            skp["explicit_kpoints_labels"][start],
            skp["explicit_kpoints_rel"][start][0],
            skp["explicit_kpoints_rel"][start][1],
            skp["explicit_kpoints_rel"][start][2],
            skp["explicit_kpoints_labels"][final],
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2]))
    #
    # Number of electrons
    #
    nelec = 0.0
    for iat in atom:
        nelec += valence_dict[iat]
    #
    # Shell scripts
    #
    structure2 = pymatgen.Structure(skp["primitive_lattice"],
                                    skp["primitive_types"], skp["primitive_positions"])
    spg_analysis = SpacegroupAnalyzer(structure2)
    middle = spg_analysis.get_ir_reciprocal_mesh(mesh=(nq[0]*2, nq[1]*2, nq[2]*2), is_shift=(0, 0, 0))
    dense = spg_analysis.get_ir_reciprocal_mesh(mesh=(nq[0]*4, nq[1]*4, nq[2]*4), is_shift=(0, 0, 0))
    print("Number of irreducible k : ", len(middle), len(dense))
    write_sh(len(middle), len(dense), len(skp["explicit_kpoints_rel"]), atom, prefix, atomwfc_dict)
    #
    # rx.in, scf.in, nscf.in, band.in , nscf_w.in, nscf_r.in
    #
    write_pwx(prefix, skp, pseudo_dir, ecutwfc, ecutrho, pseudo_dict, nq, nelec)
    #
    # ph.in, elph.in, epmat.in, phdos.in, rpa.in, scdft.in
    #
    write_ph(prefix, nq, ecutwfc, nelec)
    #
    # bands.in, pp.in, proj.in, pw2wan.in, q2r.in
    #
    write_pp(prefix)
    #
    # band.gp, {prefix}.win, respack.in, disp.in
    #
    write_wannier(prefix, skp, nelec, nq)
    #
    # openmx.in : Input file for openmx
    #
    if not os.path.isfile("openmx.in"):
        write_openmx(prefix, skp, nq)
