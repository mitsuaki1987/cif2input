from openmx import omx_pao_dict, omx_pot_dict, omx_radius_dict, omx_valence_dict
import pymatgen
import numpy
from pymatgen.core.periodic_table import get_el_sp


def write_openmx(skp, nq, rel):
    #
    # Lattice information
    #
    avec = skp["primitive_lattice"]
    pos = skp["primitive_positions"]
    nat = len(skp["primitive_types"])
    atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    typ = set(atom)
    ntyp = len(typ)
    #
    with open("openmx.dat", 'w') as f:
        print("#", file=f)
        print("# File Name", file=f)
        print("#", file=f)
        print("System.CurrentDirectory    ./", file=f)
        print("System.Name          openmx", file=f)
        print("level.of.stdout      1 #1-3", file=f)
        print("level.of.fileout     0 #0-2", file=f)
        print("data.path     ../DFT_DATA19/", file=f)
        print("HS.fileout     off   # on|off", file=f)
        print("scf.restart     off", file=f)
        #
        print("#", file=f)
        print("# Atomic Structure", file=f)
        print("#", file=f)
        print("Species.Number  %d" % (ntyp*2), file=f)
        print("<Definition.of.Atomic.Species", file=f)
        for ityp in typ:
            print(" %s  %s%s-%s  %s  %f" % (
                ityp, ityp, omx_radius_dict[str(ityp)], omx_pao_dict[str(ityp)], omx_pot_dict[str(ityp)],
                pymatgen.Element(ityp).atomic_mass), file=f)
            print("proj%s  %s%s-s1p1d1  %s" % (
                ityp, ityp, omx_radius_dict[str(ityp)], omx_pot_dict[str(ityp)]), file=f)
        print("Definition.of.Atomic.Species>", file=f)
        print("Atoms.Number  %d" % nat, file=f)
        print("Atoms.SpeciesAndCoordinates.Unit   Ang", file=f)
        print("<Atoms.SpeciesAndCoordinates", file=f)
        for iat in range(nat):
            pos2 = numpy.dot(pos[iat, :], avec)
            print("%d %s %f %f %f %f %f" % (
                iat + 1, atom[iat], pos2[0], pos2[1], pos2[2],
                omx_valence_dict[atom[iat]] * 0.5, omx_valence_dict[atom[iat]] * 0.5), file=f)
        print("Atoms.SpeciesAndCoordinates>", file=f)
        print("Atoms.UnitVectors.Unit  Ang", file=f)
        print("<Atoms.UnitVectors", file=f)
        for ii in range(3):
            print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
        print("Atoms.UnitVectors>", file=f)
        #
        print("#", file=f)
        print("# SCF or Electronic System", file=f)
        print("#", file=f)
        print("scf.XcType               GGA-PBE", file=f)
        if rel:
            print("scf.SpinPolarization        NC   # On|Off|NC", file=f)
            print("scf.SpinOrbit.Coupling      On   # On|Off", file=f)
        else:
            print("scf.SpinPolarization        Off   # On|Off|NC", file=f)
            print("scf.SpinOrbit.Coupling      off   # On|Off", file=f)
        print("<scf.SO.factor", file=f)
        for ityp in typ:
            print(" %s  s 1.0 p 1.0 d 1.0 f 1.0" % ityp, file=f)
            print("proj%s  s 1.0 p 1.0 d 1.0 f 1.0" % ityp, file=f)
        print("scf.SO.factor>", file=f)
        print("scf.ElectronicTemperature   5000", file=f)
        print("scf.maxIter                  40", file=f)
        print("scf.energycutoff           300", file=f)
        print("scf.EigenvalueSolver       band        # DC|GDC|Cluster|Band|NEGF", file=f)
        print("scf.Kgrid              %d %d %d" % (nq[0] * 2, nq[1] * 2, nq[2] * 2), file=f)
        print("scf.Mixing.Type           rmm-diisk   #Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk", file=f)
        print("scf.Init.Mixing.Weight     0.30", file=f)
        print("scf.Min.Mixing.Weight      0.001", file=f)
        print("scf.Max.Mixing.Weight      0.4", file=f)
        print("scf.Mixing.History          50", file=f)
        print("scf.Mixing.StartPulay       20", file=f)
        print("scf.Mixing.EveryPulay       5", file=f)
        print("#scf.Kerker.factor  ", file=f)
        print("scf.criterion       %e" % (float(nat)*1.0e-10), file=f)
        print("scf.partialCoreCorrection    on", file=f)
        print("scf.partialCoreCorrection    on", file=f)
        print("scf.system.charge            0.0", file=f)
        #
        print("#", file=f)
        print("# DFT+U", file=f)
        print("#", file=f)
        print("scf.Hubbard.U          off", file=f)
        print("scf.Hubbard.Occupation     dual   # onsite|full|dual", file=f)
        print("<Hubbard.U.values", file=f)
        for ityp in typ:
            print(" %s  1s 0.0 1p 0.0 1d 0.0" % ityp, file=f)
            print("proj%s  1s 0.0 1p 0.0 1d 0.0" % ityp, file=f)
        print("Hubbard.U.values>", file=f)
        #
        print("#", file=f)
        print("# 1D FFT", file=f)
        print("#", file=f)
        print("1DFFT.EnergyCutoff       3600", file=f)
        print("1DFFT.NumGridK            900", file=f)
        print("1DFFT.NumGridR            900", file=f)
        #
        print("#", file=f)
        print("# van der Waals", file=f)
        print("#", file=f)
        print("scf.DFTD            off", file=f)
        print("version.dftD         2  # 2|3", file=f)
        print("DFTD.scale6            1   # default=0.75|1.0 (for DFT-D2|DFT-D3)", file=f)
        print("DFTD.scale8       0.7875   # default=0.722|0.7875 (for PBE with zero|bj damping)", file=f)
        print("DFTD.sr6           1.217", file=f)
        print("DFTD.a1           0.4289", file=f)
        print("DFTD.a2           4.4407", file=f)
        print("DFTD.cncut_dftD              40            # default=40 (DFTD.Unit)", file=f)
        print("DFTD.Unit                   Ang", file=f)
        print("DFTD.rcut_dftD             100.0", file=f)
        print("DFTD.d                      20.0", file=f)
        print("DFTD.scale6                 0.75", file=f)
        print("DFTD.IntDirection          1 1 1", file=f)
        print("<DFTD.periodicity", file=f)
        for iat in range(nat):
            print("%d 1" % (iat + 1), file=f)
        print("DFTD.periodicity>", file=f)
        #
        print("#", file=f)
        print("# MD and Structure optimization", file=f)
        print("#", file=f)
        print("orbitalOpt.Force.Skip       off", file=f)
        print("scf.stress.tensor  off  # on|off, default=off", file=f)
        print("MD.Type    NOMD", file=f)
        print("# NOMD|Opt|NVE|NVT_VS|NVT_VS2|NVT_NH", file=f)
        print("# Opt|DIIS|BFGS|RF|EF|", file=f)
        print("# OptC1|OptC2|OptC3|OptC4|OptC5|RFC5", file=f)
        print("# NEB", file=f)
        print("MD.Opt.DIIS.History     3", file=f)
        print("MD.Opt.StartDIIS      30", file=f)
        print("MD.Opt.EveryDIIS       200", file=f)
        print("MD.maxIter             100", file=f)
        print("MD.Opt.criterion      0.0005", file=f)
        print("MD.Opt.Init.Hessian        Schlegel  # Schlegel|iden", file=f)
        print("<MD.Fixed.XYZ", file=f)
        for iat in range(nat):
            print("%d 0 0 0" % (iat + 1), file=f)
        print("MD.Fixed.XYZ>", file=f)
        print("<MD.Fixed.Cell.Vectors", file=f)
        print("0 0 0", file=f)
        print("0 0 0", file=f)
        print("0 0 0", file=f)
        print("MD.Fixed.Cell.Vectors>", file=f)
        print("<MD.TempControl", file=f)
        print("3", file=f)
        print("100   2  1000.0  0.0  ", file=f)
        print("400  10   700.0  0.4  ", file=f)
        print("700  40   500.0  0.7  ", file=f)
        print("MD.TempControl>", file=f)
        print("<MD.Init.Velocity", file=f)
        for iat in range(nat):
            print("%d 0.0 0.0 0.0" % (iat + 1), file=f)
        print("MD.Init.Velocity>", file=f)
        #
        print("#", file=f)
        print("# Band dispersion", file=f)
        print("#", file=f)
        print("Band.dispersion              off", file=f)
        print("Band.Nkpath  %d" % len(skp["path"]), file=f)
        print("<Band.kpath", file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            final = skp["explicit_segments"][ipath][1] - 1
            print("%d  %f %f %f %f %f %f %s %s" % (
                final - start + 1,
                skp["explicit_kpoints_rel"][start][0],
                skp["explicit_kpoints_rel"][start][1],
                skp["explicit_kpoints_rel"][start][2],
                skp["explicit_kpoints_rel"][final][0],
                skp["explicit_kpoints_rel"][final][1],
                skp["explicit_kpoints_rel"][final][2],
                skp["explicit_kpoints_labels"][start],
                skp["explicit_kpoints_labels"][final]),
                  file=f)
        print("Band.kpath>", file=f)
        #
        print("#", file=f)
        print("# Wannier", file=f)
        print("#", file=f)
        print("Wannier.Func.Calc off", file=f)
        print("Wannier.Func.Num %d" % (nat*9), file=f)
        print("Wannier.Outer.Window.Bottom  -1.5", file=f)
        print("Wannier.Outer.Window.Top      7.0", file=f)
        print("Wannier.Inner.Window.Bottom  -1.5", file=f)
        print("Wannier.Inner.Window.Top      1.2", file=f)
        print("Wannier.Initial.Projectors.Unit ANG", file=f)
        print("<Wannier.Initial.Projectors", file=f)
        for iat in range(nat):
            pos2 = numpy.dot(pos[iat, :], avec)
            for prj in "s", "px", "py", "pz", "dxy", "dyz", "dxz", \
                       "dx2-y2", "dz2":
                print("proj%s-%s %f %f %f  0.0 0.0 1.0  1.0 0.0 0.0" % (
                      atom[iat], prj,
                      pos2[0], pos2[1], pos2[2]), file=f)
        print("Wannier.Initial.Projectors>", file=f)
        print("Wannier.Interpolated.Bands             off", file=f)
        print("Wannier.Function.Plot                  off", file=f)
        print("Wannier.Function.Plot.SuperCells      1 1 1", file=f)
        print("Wannier.Minimizing.Scheme             0", file=f)
        print("Wannier.Minimizing.StepLength        2.0", file=f)
        print("Wannier.Minimizing.Secant.Steps       5", file=f)
        print("Wannier.Minimizing.Secant.StepLength 2.0", file=f)
        print("Wannier.Minimizing.Conv.Criterion   1e-8", file=f)
        print("Wannier.Minimizing.Max.Steps         200", file=f)
        print("Wannier.Readin.Overlap.Matrix       on", file=f)
        print("Wannier.Dis.SCF.Max.Steps         200", file=f)
        print("Wannier.Dis.Conv.Criterion        1e-8", file=f)
        print("Wannier.Dis.Mixing.Para           0.5", file=f)
        print("Wannier.MaxShells          12", file=f)
        print("Wannier.Kgrid     %d %d %d" % (nq[0], nq[1], nq[2]), file=f)
        #
        print("#", file=f)
        print("# DOS", file=f)
        print("#", file=f)
        print("Dos.fileout      off", file=f)
        print("Dos.Erange       -20.0  20.0", file=f)
        print("Dos.Kgrid     %d %d %d" % (nq[0] * 4, nq[1] * 4, nq[2] * 4), file=f)
        print("FermiSurfer.fileout         on", file=f)
        print("OpticalConductivity.fileout    off", file=f)
        #
        print("#", file=f)
        print("# Energy Decomposition", file=f)
        print("#", file=f)
        print("Energy.Decomposition      off        # on|off", file=f)
        #
        print("#", file=f)
        print("# Orbital optimization", file=f)
        print("#", file=f)
        print("orbitalOpt.Method      off   # Off|Species|Atoms", file=f)
        print("orbitalOpt.Opt.Method     EF     # DIIS|EF", file=f)
        print("orbitalOpt.SD.step       0.001", file=f)
        print("orbitalOpt.HistoryPulay     15", file=f)
        print("orbitalOpt.StartPulay       1", file=f)
        print("orbitalOpt.scf.maxIter     40", file=f)
        print("orbitalOpt.Opt.maxIter    100", file=f)
        print("orbitalOpt.per.MDIter     1000000", file=f)
        print("orbitalOpt.criterion      1.0e-4", file=f)
        print("CntOrb.fileout       off        # on|off", file=f)
        print("Num.CntOrb.Atom    %d" % nat, file=f)
        print("<Atoms.Cont.Orbitals", file=f)
        for iat in range(nat):
            print("%d" % iat, file=f)
        print("Atoms.Cont.Orbitals>", file=f)
        #
        print("#", file=f)
        print("# Order-N", file=f)
        print("#", file=f)
        print("orderN.HoppingRanges    6.0", file=f)
        print("orderN.KrylovH.order    400", file=f)
        print("orderN.Exact.Inverse.S   on    #on|off", file=f)
        print("orderN.KrylovS.order   1600", file=f)
        print("orderN.Recalc.Buffer on #on|off", file=f)
        print("orderN.Expand.Core  on  #on|off", file=f)
        #
        print("#", file=f)
        print("# Electric Field", file=f)
        print("#", file=f)
        print("scf.Electric.Field   0.0 0.0 0.0", file=f)
        #
        print("#", file=f)
        print("# Natural population analysis ", file=f)
        print("#", file=f)
        print("NBO.switch   off # on1|on2", file=f)
        print("NBO.Num.CenterAtoms     5", file=f)
        print("<NBO.CenterAtoms", file=f)
        print("269", file=f)
        print("304", file=f)
        print("323", file=f)
        print("541", file=f)
        print("574", file=f)
        print("NBO.CenterAtoms>", file=f)
        #
        print("#", file=f)
        print("# Magnetic field", file=f)
        print("#", file=f)
        print("scf.Constraint.NC.Spin       off    # on|on2|off", file=f)
        print("scf.Constraint.NC.Spin.v    0.0", file=f)
        print("scf.NC.Zeeman.Spin       off        # on|off", file=f)
        print("scf.NC.Mag.Field.Spin    0.0", file=f)
        print("scf.NC.Zeeman.Orbital      off        # on|off", file=f)
        print("scf.NC.Mag.Field.Orbital   0.0", file=f)
        #
        print("#", file=f)
        print("# ESM", file=f)
        print("#", file=f)
        print("ESM.switch    off    # off, on1=v|v|v, on2=m|v|m, on3=v|v|m, on4=on2+EF", file=f)
        print("ESM.buffer.range      10.0", file=f)
        print("ESM.potential.diff     0.0", file=f)
        print("ESM.wall.position        10.0", file=f)
        print("ESM.wall.height        100.0", file=f)
        #
        print("#", file=f)
        print("# NEB", file=f)
        print("#", file=f)
        print("MD.NEB.Number.Images     10", file=f)
        print("MD.NEB.Spring.Const      0.1", file=f)
        print("<NEB.Atoms.SpeciesAndCoordinates", file=f)
        for iat in range(nat):
            pos2 = numpy.dot(pos[iat, :], avec)
            print("%d %s %f %f %f %f %f" % (
                iat + 1, atom[iat], pos2[0], pos2[1], pos2[2],
                omx_valence_dict[atom[iat]] * 0.5, omx_valence_dict[atom[iat]] * 0.5), file=f)
        print("NEB.Atoms.SpeciesAndCoordinates>", file=f)
        #
        print("#", file=f)
        print("# STM Image", file=f)
        print("#", file=f)
        print("partial.charge     off    # on|off", file=f)
        print("partial.charge.energy.window   0.0", file=f)
        #
        print("#", file=f)
        print("# Band Unfolding", file=f)
        print("#", file=f)
        print("Unfolding.Electronic.Band      off    # on|off", file=f)
        print("Unfolding.LowerBound        -10", file=f)
        print("Unfolding.UpperBound          10", file=f)
        print("<Unfolding.kpoint", file=f)
        final = 0
        ii = 1
        print("%s %f %f %f" % (
            skp["explicit_kpoints_labels"][final],
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2]),
              file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start == final:
                final = skp["explicit_segments"][ipath][1] - 1
                ii += 1
                print("%s %f %f %f" % (
                    skp["explicit_kpoints_labels"][final],
                    skp["explicit_kpoints_rel"][final][0],
                    skp["explicit_kpoints_rel"][final][1],
                    skp["explicit_kpoints_rel"][final][2]),
                      file=f)
            else:
                break
        print("Unfolding.kpoint>", file=f)
        print("Unfolding.Nkpoint     %d" % ii, file=f)
        print("Unfolding.desired_totalnkpt    100", file=f)
        #
        print("#", file=f)
        print("# Draw Kohn-Sham orbital", file=f)
        print("#", file=f)
        print("MO.fileout off", file=f)
        print("num.HOMOs  1", file=f)
        print("num.LUMOs  1", file=f)
        print("MO.Nkpoint  1", file=f)
        print("<MO.kpoint", file=f)
        print("0.0  0.0  0.0", file=f)
        print("MO.kpoint>", file=f)
        #
        print("#", file=f)
        print("# NEGF", file=f)
        print("#", file=f)
        print("NEGF.output_hks    off", file=f)
        print("NEGF.filename.hks  left.hks", file=f)
        #
        print("NEGF.filename.hks.l   left.hks", file=f)
        print("NEGF.filename.hks.r   right.hks", file=f)
        print("LeftLeadAtoms.Number  %d" % nat, file=f)
        print("<LeftLeadAtoms.SpeciesAndCoordinates         ", file=f)
        for iat in range(nat):
            pos2 = numpy.dot(pos[iat, :], avec)
            print("%d %s %f %f %f %f %f" % (
                iat + 1, atom[iat], pos2[0], pos2[1], pos2[2],
                omx_valence_dict[atom[iat]] * 0.5, omx_valence_dict[atom[iat]] * 0.5), file=f)
        print("LeftLeadAtoms.SpeciesAndCoordinates>", file=f)
        print("RightLeadAtoms.Number  %d" % nat, file=f)
        print("<RightLeadAtoms.SpeciesAndCoordinates", file=f)
        for iat in range(nat):
            pos2 = numpy.dot(pos[iat, :], avec)
            print("%d %s %f %f %f %f %f" % (
                iat + 1, atom[iat], pos2[0], pos2[1], pos2[2],
                omx_valence_dict[atom[iat]] * 0.5, omx_valence_dict[atom[iat]] * 0.5), file=f)
        print("RightLeadAtoms.SpeciesAndCoordinates>", file=f)
        #
        print("NEGF.Num.Poles             150", file=f)
        print("NEGF.scf.Kgrid             1 1", file=f)
        print("NEGF.bias.voltage          0.0", file=f)
        print("NEGF.bias.neq.im.energy    0.01", file=f)
        print("NEGF.bias.neq.energy.step  0.02", file=f)
        print("NEGF.scf.Iter.Band          6", file=f)
        print("NEGF.Poisson.Solver       FD     # FD|FFT", file=f)
        print("NEGF.gate.voltage   0.0", file=f)
        #
        print("NEGF.tran.SCF.skip  off", file=f)
        print("NEGF.tran.Analysis         on", file=f)
        print("NEGF.tran.CurrentDensity   on", file=f)
        print("NEGF.tran.Channel          on", file=f)
        print("NEGF.tran.energyrange -10 10 1.0e-3", file=f)
        print('NEGF.tran.energydiv        200', file=f)
        print("NEGF.tran.Kgrid            1 1", file=f)
        print("NEGF.Channel.Nkpoint        1", file=f)
        print("<NEGF.Channel.kpoint", file=f)
        print("0.0  0.0", file=f)
        print("NEGF.Channel.kpoint>", file=f)
        print("NEGF.Channel.Nenergy        1", file=f)
        print("<NEGF.Channel.energy", file=f)
        print("0.0", file=f)
        print("NEGF.Channel.energy>", file=f)
        print("NEGF.Channel.Num    5", file=f)
        print("NEGF.Dos.energyrange     -10.0 10.0 5.0e-3", file=f)
        print("NEGF.Dos.energy.div        200", file=f)
        print("NEGF.Dos.Kgrid             1 1", file=f)
        #
        print("NEGF.tran.interpolate         off    # on|off", file=f)
        print("NEGF.tran.interpolate.file1  start.tranb", file=f)
        print("NEGF.tran.interpolate.file2  end.tranb", file=f)
        print("NEGF.tran.interpolate.coes   1.0 0.0", file=f)
