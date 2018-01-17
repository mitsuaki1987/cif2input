#!/usr/bin/python3
import sys
import numpy
import os

args = sys.argv
nfile = len(args) - 1

for ifile in range(nfile):
    print(args[ifile+1], " -> ", args[ifile+1]+"2.xsf")

    f = open(args[ifile+1], 'r')
    #
    # Skip
    #
    for ii in range(6):
        line = f.readline()
    #
    # Direct lattice vector
    #
    avec = numpy.zeros((3, 3), numpy.float_)
    for ii in range(3):
        line = f.readline()
        sp = line.split()
        for jj in range(3):
            avec[ii, jj] = float(sp[jj])
    #
    # Skip
    #
    for ii in range(5):
        line = f.readline()
    #
    # Atomic position
    #
    line = f.readline()
    sp = line.split()
    natom = int(sp[0])
    pos = numpy.zeros((natom, 3), numpy.float_)
    label = [""]*natom
    for iatom in range(natom):
        line = f.readline()
        sp = line.split()
        label[iatom] = sp[0]
        for ii in range(3):
            pos[iatom, ii] = float(sp[ii+1])
    #
    # Skip
    #
    for ii in range(5):
        line = f.readline()
    #
    # grid and origin
    #
    line = f.readline()
    sp = line.split()
    ng = numpy.zeros(3, numpy.int_)
    for ii in range(3):
        ng[ii] = int(sp[ii])
    line = f.readline()
    sp = line.split()
    origin = numpy.zeros(3, numpy.float_)
    for ii in range(3):
        origin[ii] = float(sp[ii])
    for ii in range(3):
        origin[:] += avec[ii, :]
    f.close()
    #
    # Print file
    #
    f = open(args[ifile+1]+"2.xsf", 'w')
    print("CRYSTAL", file=f)
    print("PRIMVEC", file=f)
    for ii in range(3):
        print("%f %f %f" % (3.0*avec[ii, 0], 3.0*avec[ii, 1], 3.0*avec[ii, 2]), file=f)
    print("PRIMCOORD", file=f)
    print("%d 1" % (natom*27), file=f)
    pos2 = numpy.zeros(3, numpy.float_)
    for i0 in range(3):
        for i1 in range(3):
            for i2 in range(3):
                for iatom in range(natom):
                    pos2[:] = pos[iatom, :]+i0*avec[0, :]+i1*avec[1, :]+i2*avec[2, :]
                    print("%s %f %f %f" % (label[iatom], pos2[0], pos2[1], pos2[2]), file=f)
    print("BEGIN_BLOCK_DATAGRID_3D", file=f)
    print("3D_field", file=f)
    print("BEGIN_DATAGRID_3D_UNKNOWN", file=f)
    print("%d %d %d" % (ng[0], ng[1], ng[2]), file=f)
    print("%f %f %f" % (origin[0], origin[1], origin[2]), file=f)
    f.close()
    #
    #
    #
    check = os.system('sed -e \"1,{0}d\" {1} >> {2}'.format(22+natom, args[ifile+1], args[ifile+1]+"2.xsf"))
