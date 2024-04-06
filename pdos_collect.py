#!/usr/bin/python3
from pymatgen.core.structure import Structure
import numpy
import sys
import spglib
import libtetrabz


def main():
    #
    args = sys.argv
    with open(str(args[1]), "r") as f:
        input_list = f.readlines()
    #
    # Read previous result
    #
    for input_file in input_list:
        #
        input_file = input_file.strip("\n")
        prefix = input_file.split("/")[-1].split(".")[0]
        #
        # Crystal structure and symmetry
        #
        crystal = Structure.from_file("ml/"+prefix + ".xsf")
        avec = crystal.lattice.matrix
        pos = numpy.array([atom.frac_coords for atom in crystal])
        spc = [atom.specie.number for atom in crystal]
        rot = spglib.get_symmetry((avec, pos, spc))["rotations"]
        eq_atom = spglib.get_symmetry((avec, pos, spc))["equivalent_atoms"]
        bvec = crystal.lattice.reciprocal_lattice.matrix
        nat = len(crystal)
        #
        # Fermi energy
        #
        alat = 0.0
        efermi = 0.0
        with open("nscfout/nscf_" + prefix + ".out") as f:
            filelines = f.readlines()
        for fileline in filelines:
            if fileline[0:32] == "     lattice parameter (alat)  =":
                alat = float(fileline[32:45]) * 0.529177249
            elif fileline[0:24] == "     the Fermi energy is":
                efermi = float(fileline[24:35])
        #
        # k-point grid
        #
        nk = numpy.zeros([3], dtype=numpy.int_)
        with open("nscfin/nscf_" + prefix + ".in") as f:
            filelines = f.readlines()
        for iline in range(len(filelines)):
            if filelines[iline][0:18] == "K_POINTS automatic":
                nkstr = filelines[iline + 1].split(" ")
                nk = numpy.array([nkstr[1], nkstr[2], nkstr[3]], dtype=numpy.int_)
        #
        # k-point in irreducible BZ (IBZ), band energy, projection in IBZ
        #
        nbnd = 0
        nk_ibz = 0
        wfc2atom = numpy.empty([1], dtype=numpy.int_)
        with open("projout/pdos_" + prefix + ".out") as f:
            filelines = f.readlines()
        nline = len(filelines)
        for iline in range(nline):
            if filelines[iline] == "  Problem Sizes \n":
                natomwfc = int(filelines[iline + 1].split("=")[1])
                nbnd = int(filelines[iline + 2].split("=")[1])
                nk_ibz = int(int(filelines[iline + 3].split("=")[1]) / 2)
                wfc2atom = numpy.empty([natomwfc], dtype=numpy.int_)
                for iwfc in range(natomwfc):
                    wfc2atom[iwfc] = int(filelines[iline + 11 + iwfc][22:26]) - 1
                break

        kvec = numpy.empty([nk_ibz, 3], dtype=numpy.int_)
        # band energy and projection in IBZ
        eig_ibz = numpy.empty([nk_ibz, nbnd * 2], dtype=numpy.float_)
        proj_ibz = numpy.zeros([nk_ibz, nbnd * 2, nat], dtype=numpy.float_)
        #
        ik = 0
        ispin = 0
        for iline in range(nline):
            if filelines[iline][0:4] == " k =":
                kvec0 = numpy.array([filelines[iline][5:19],
                                     filelines[iline][19:33],
                                     filelines[iline][33:47]], dtype=numpy.float_)
                kvec0 = numpy.dot(avec, kvec0) / alat
                kvec0 = kvec0 * nk
                ikvec = [round(kvec1) for kvec1 in kvec0]
                kvec[ik, :] = numpy.array(ikvec)

                ibnd = 0
                for jline in range(iline + 1, nline):
                    if filelines[jline][0:7] == "==== e(":
                        eig0 = float(filelines[jline][14:26])
                        eig_ibz[ik, ibnd + ispin * nbnd] = eig0

                        proj_str = ""
                        for kline in range(jline + 1, nline):
                            if filelines[kline][0:13] == "    |psi|^2 =":
                                break
                            else:
                                proj_str += filelines[kline][10:].strip("\n")
                        proj_str = proj_str.replace("]", "")
                        proj_str = proj_str.split("+")
                        proj_str = [proj_str_wfc.split("*[#") for proj_str_wfc in proj_str]
                        for proj_str_wfc in proj_str:
                            if len(proj_str_wfc) != 2:
                                continue
                            iwfc = int(proj_str_wfc[1]) - 1
                            proj_wfc = float(proj_str_wfc[0])
                            proj_ibz[ik, ibnd + ispin * nbnd, wfc2atom[iwfc]] += proj_wfc

                        ibnd += 1
                        if ibnd == nbnd:
                            break

                ik += 1
                if ik == nk_ibz:
                    ik = 0
                    ispin = 1
        #
        # Projection becomes the same for equivalent atoms
        #
        proj2 = numpy.zeros([nk_ibz, 2 * nbnd], dtype=numpy.float_)
        for iat in range(nat):
            proj2[:, :] = 0.0
            nat0 = 0
            for jat in range(nat):
                if eq_atom[jat] == iat:
                    nat0 += 1
                    proj2[:, :] += proj_ibz[:, :, jat]
            if nat0 == 0:
                proj_ibz[:, :, iat] = proj_ibz[:, :, eq_atom[iat]]
            else:
                proj_ibz[:, :, iat] = proj2[:, :] / nat0
        #
        # If insulator, shift Fermi level to the valence top
        #
        occ = numpy.empty([2 * nbnd], dtype=numpy.int_)
        for ibnd in range(2 * nbnd):
            emax = max(eig_ibz[:, ibnd])
            emin = min(eig_ibz[:, ibnd])
            if emax < efermi:
                occ[ibnd] = 1
            elif efermi < emin:
                occ[ibnd] = -1
            else:
                occ[ibnd] = 0
        if not numpy.any(occ == 0):
            valence_top = -9999.0
            conduct_btm = 9999.0
            for ibnd in range(2 * nbnd):
                if occ[ibnd] == 1:
                    valence_top = max(valence_top, max(eig_ibz[:, ibnd]))
                elif occ[ibnd] == -1:
                    conduct_btm = min(conduct_btm, min(eig_ibz[:, ibnd]))
            efermi = 0.5 * (valence_top + conduct_btm)
        #
        # energy and projection : IBZ -> full BZ with symmetry operator
        #
        eig_fbz = numpy.empty([nk[0], nk[1], nk[2], 2 * nbnd], dtype=numpy.float_)
        proj_fbz = numpy.empty([nk[0], nk[1], nk[2], 2 * nbnd, nat], dtype=numpy.float_)
        kdone = numpy.zeros([nk[0], nk[1], nk[2]], dtype=numpy.int_)
        for ik in range(nk_ibz):
            for rot0 in rot:
                ikv1 = numpy.dot(kvec[ik, :], rot0)
                ikv1 = ikv1 % nk
                kdone[ikv1[0], ikv1[1], ikv1[2]] += 1
                eig_fbz[ikv1[0], ikv1[1], ikv1[2], :] = eig_ibz[ik, :]
                proj_fbz[ikv1[0], ikv1[1], ikv1[2], :, :] = proj_ibz[ik, :, :]
                #
                # time-reversal symmetry
                #
                ikv1 = - ikv1
                kdone[ikv1[0], ikv1[1], ikv1[2]] += 1
                eig_fbz[ikv1[0], ikv1[1], ikv1[2], :] = eig_ibz[ik, :]
                proj_fbz[ikv1[0], ikv1[1], ikv1[2], :, :] = proj_ibz[ik, :, :]
        #
        # Weights for partial density of states at E_F and averaged PDOS
        # computed with optimized tetrahedron method
        #
        de = numpy.array([0.0, 0.01, 0.03, 0.1, 0.3, 1.0, -1.0, -0.3, -0.1, -0.03, -0.01])
        e0 = de + efermi
        wght = libtetrabz.intdos(bvec, eig_fbz[:, :, :, :], e0)
        wght0 = libtetrabz.dos(bvec, eig_fbz[:, :, :, :], numpy.array([efermi]))
        for ie in range(1, 6):
            wght[:, :, :, :, ie] = (wght[:, :, :, :, ie] - wght[:, :, :, :, -ie]) / (de[ie] - de[-ie])
        wght[:, :, :, :, 0] = wght0[:, :, :, :, 0]
        #
        # Total DOS
        #
        for ie in range(6):
            with open(str(args[1])+"_dos"+str(de[ie])+".csv", "a") as f:
                print(prefix + "," + str(numpy.sum(wght[:, :, :, :, ie])), file=f)
        #
        # Partial DOS
        #
        for ie in range(6):
            with open(str(args[1]) + "_pdos" + str(de[ie]) + ".csv", "a") as f:
                print(prefix, end="", file=f)
                for iat in range(nat):
                    print("," + str(numpy.sum(wght[:, :, :, :, ie] * proj_fbz[:, :, :, :, iat])), end="", file=f)
                print("", file=f)


main()
