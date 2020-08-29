#!/usr/bin/python3
import sys
import pymatgen
import seekpath
import numpy

args = sys.argv
structure = pymatgen.Structure.from_file(args[1])
structure.remove_oxidation_states()
skp = seekpath.get_explicit_k_path((structure.lattice.matrix, numpy.array(structure.frac_coords),
                                   [pymatgen.Element(str(spc)).number for spc in structure.species]),
                                   reference_distance=0.1)

print("Band.dispersion     on")
print("<Band.KPath.UnitCell")
for ii in range(3):
    print(" %f %f %f" % (skp["primitive_lattice"][ii, 0],
                         skp["primitive_lattice"][ii, 1],
                         skp["primitive_lattice"][ii, 2]))
print("Band.KPath.UnitCell>")
print("Band.Nkpath  %d" % len(skp["path"]))
print("<Band.kpath")
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
          skp["explicit_kpoints_labels"][final]))
print("Band.kpath>")

