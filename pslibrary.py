import numpy


pseudo_dict = {
    "Ac": "Ac.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Ag": "Ag.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Al": "Al.pbe-n-rrkjus_psl.0.1.UPF",
    "Ar": "Ar.pbe-n-rrkjus_psl.1.0.0.UPF",
    "As": "As.pbe-n-rrkjus_psl.0.2.UPF",
    "At": "At.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Au": "Au.pbe-dn-rrkjus_psl.0.1.UPF",
    "B": "B.pbe-n-rrkjus_psl.0.1.UPF",
    "Ba": "Ba.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Be": "Be.pbe-n-rrkjus_psl.0.2.UPF",
    "Bi": "Bi.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Br": "Br.pbe-n-rrkjus_psl.0.2.UPF",
    "C": "C.pbe-n-rrkjus_psl.0.1.UPF",
    "Ca": "Ca.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Cd": "Cd.pbe-dn-rrkjus_psl.0.2.UPF",
    "Ce": "Ce.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Cl": "Cl.pbe-n-rrkjus_psl.0.1.UPF",
    "Co": "Co.pbe-spn-rrkjus_psl.0.3.1.UPF",
    "Cr": "Cr.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Cs": "Cs.pbe-spn-rrkjus_psl.0.2.3.UPF",
    "Cu": "Cu.pbe-dn-rrkjus_psl.0.2.UPF",
    "Dy": "Dy.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Er": "Er.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Eu": "Eu.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "F": "F.pbe-n-rrkjus_psl.0.1.UPF",
    "Fe": "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF",
    "Fr": "Fr.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Ga": "Ga.pbe-dn-rrkjus_psl.0.2.UPF",
    "Gd": "Gd.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Ge": "Ge.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "H": "H.pbe-rrkjus_psl.0.1.UPF",
    "He": "He.pbe-rrkjus_psl.1.0.0.UPF",
    "Hf": "Hf.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Hg": "Hg.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Ho": "Ho.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "I": "I.pbe-n-rrkjus_psl.0.2.UPF",
    "In": "In.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Ir": "Ir.pbe-n-rrkjus_psl.0.2.3.UPF",
    "K": "K.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Kr": "Kr.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "La": "La.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Li": "Li.pbe-s-rrkjus_psl.0.2.1.UPF",
    "Lu": "Lu.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Mg": "Mg.pbe-n-rrkjus_psl.0.3.0.UPF",
    "Mn": "Mn.pbe-spn-rrkjus_psl.0.3.1.UPF",
    "Mo": "Mo.pbe-spn-rrkjus_psl.0.2.UPF",
    "N": "N.pbe-n-rrkjus_psl.0.1.UPF",
    "Na": "Na.pbe-spn-rrkjus_psl.0.2.UPF",
    "Nb": "Nb.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Nd": "Nd.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Ne": "Ne.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Ni": "Ni.pbe-n-rrkjus_psl.0.1.UPF",
    "Np": "Np.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "O": "O.pbe-n-rrkjus_psl.0.1.UPF",
    "Os": "Os.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "P": "P.pbe-n-rrkjus_psl.0.1.UPF",
    "Pa": "Pa.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Pb": "Pb.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Pd": "Pd.pbe-n-rrkjus_psl.0.2.2.UPF",
    "Pm": "Pm.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Po": "Po.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Pr": "Pr.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Pt": "Pt.pbe-n-rrkjus_psl.0.1.UPF",
    "Pu": "Pu.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Ra": "Ra.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Rb": "Rb.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Re": "Re.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Rh": "Rh.pbe-spn-rrkjus_psl.0.2.3.UPF",
    "Rn": "Rn.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Ru": "Ru.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "S": "S.pbe-n-rrkjus_psl.0.1.UPF",
    "Sb": "Sb.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Sc": "Sc.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Se": "Se.pbe-n-rrkjus_psl.0.2.UPF",
    "Si": "Si.pbe-n-rrkjus_psl.0.1.UPF",
    "Sm": "Sm.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Sn": "Sn.pbe-dn-rrkjus_psl.0.2.UPF",
    "Sr": "Sr.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Ta": "Ta.pbe-spn-rrkjus_psl.0.2.UPF",
    "Tb": "Tb.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Tc": "Tc.pbe-spn-rrkjus_psl.0.3.0.UPF",
    "Te": "Te.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Th": "Th.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Ti": "Ti.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Tl": "Tl.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Tm": "Tm.pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "U": "U.pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "V": "V.pbe-spnl-rrkjus_psl.1.0.0.UPF",
    "W": "W.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Xe": "Xe.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Y": "Y.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Yb": "Yb.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Zn": "Zn.pbe-dnl-rrkjus_psl.1.0.0.UPF",
    "Zr": "Zr.pbe-spn-rrkjus_psl.1.0.0.UPF"
    }

ecutwfc_dict = {
    "Ac": 48.0,
    "Ag": 45.0,
    "Al": 36.0,
    "Ar": 48.0,
    "As": 20.0,
    "At": 50.0,
    "Au": 26.0,
    "B": 48.0,
    "Ba": 37.0,
    "Be": 30.0,
    "Bi": 45.0,
    "Br": 25.0,
    "C": 37.0,
    "Ca": 45.0,
    "Cd": 31.0,
    "Ce": 37.0,
    "Cl": 36.0,
    "Co": 60.0,
    "Cr": 48.0,
    "Cs": 30.0,
    "Cu": 33.0,
    "Dy": 36.0,
    "Er": 37.0,
    "Eu": 35.0,
    "F": 48.0,
    "Fe": 64.0,
    "Fr": 63.0,
    "Ga": 33.0,
    "Gd": 35.0,
    "Ge": 38.0,
    "H": 46.0,
    "He": 46.0,
    "Hf": 46.0,
    "Hg": 29.0,
    "Ho": 36.0,
    "I": 28.0,
    "In": 48.0,
    "Ir": 27.0,
    "K": 41.0,
    "Kr": 54.0,
    "La": 75.0,
    "Li": 78.0,
    "Lu": 42.0,
    "Mg": 13.0,
    "Mn": 46.0,
    "Mo": 48.0,
    "N": 39.0,
    "Na": 41.0,
    "Nb": 47.0,
    "Nd": 38.0,
    "Ne": 54.0,
    "Ni": 41.0,
    "Np": 111.0,
    "O": 47.0,
    "Os": 54.0,
    "P": 17.0,
    "Pa": 91.0,
    "Pb": 40.0,
    "Pd": 32.0,
    "Pm": 34.0,
    "Po": 47.0,
    "Pr": 41.0,
    "Pt": 39.0,
    "Pu": 45.0,
    "Ra": 78.0,
    "Rb": 35.0,
    "Re": 51.0,
    "Rh": 41.0,
    "Rn": 49.0,
    "Ru": 52.0,
    "S": 17.0,
    "Sb": 27.0,
    "Sc": 50.0,
    "Se": 21.0,
    "Si": 38.0,
    "Sm": 34.0,
    "Sn": 25.0,
    "Sr": 39.0,
    "Ta": 47.0,
    "Tb": 35.0,
    "Tc": 62.0,
    "Te": 34.0,
    "Th": 103.0,
    "Ti": 52.0,
    "Tl": 45.0,
    "Tm": 38.0,
    "U": 93.0,
    "V": 48.0,
    "W": 50.0,
    "Xe": 45.0,
    "Y": 38.0,
    "Yb": 39.0,
    "Zn": 44.0,
    "Zr": 46.0
}

ecutrho_dict = {
    "Ac": 328.0,
    "Ag": 181.0,
    "Al": 145.0,
    "Ar": 225.0,
    "As": 103.0,
    "At": 471.0,
    "Au": 909.0,
    "B": 216.0,
    "Ba": 176.0,
    "Be": 197.0,
    "Bi": 455.0,
    "Br": 108.0,
    "C": 318.0,
    "Ca": 274.0,
    "Cd": 156.0,
    "Ce": 179.0,
    "Cl": 146.0,
    "Co": 445.0,
    "Cr": 457.0,
    "Cs": 240.0,
    "Cu": 150.0,
    "Dy": 265.0,
    "Er": 265.0,
    "Eu": 259.0,
    "F": 316.0,
    "Fe": 782.0,
    "Fr": 506.0,
    "Ga": 180.0,
    "Gd": 227.0,
    "Ge": 237.0,
    "H": 221.0,
    "He": 203.0,
    "Hf": 264.0,
    "Hg": 119.0,
    "Ho": 268.0,
    "I": 113.0,
    "In": 190.0,
    "Ir": 361.0,
    "K": 277.0,
    "Kr": 252.0,
    "La": 537.0,
    "Li": 355.0,
    "Lu": 333.0,
    "Mg": 87.0,
    "Mn": 244.0,
    "Mo": 407.0,
    "N": 264.0,
    "Na": 163.0,
    "Nb": 266.0,
    "Nd": 202.0,
    "Ne": 264.0,
    "Ni": 236.0,
    "Np": 551.0,
    "O": 323.0,
    "Os": 268.0,
    "P": 79.0,
    "Pa": 446.0,
    "Pb": 158.0,
    "Pd": 133.0,
    "Pm": 198.0,
    "Po": 455.0,
    "Pr": 309.0,
    "Pt": 401.0,
    "Pu": 351.0,
    "Ra": 523.0,
    "Rb": 259.0,
    "Re": 274.0,
    "Rh": 440.0,
    "Rn": 213.0,
    "Ru": 353.0,
    "S": 77.0,
    "Sb": 136.0,
    "Sc": 402.0,
    "Se": 95.0,
    "Si": 151.0,
    "Sm": 195.0,
    "Sn": 126.0,
    "Sr": 262.0,
    "Ta": 387.0,
    "Tb": 229.0,
    "Tc": 832.0,
    "Te": 148.0,
    "Th": 413.0,
    "Ti": 575.0,
    "Tl": 210.0,
    "Tm": 316.0,
    "U": 683.0,
    "V": 645.0,
    "W": 475.0,
    "Xe": 212.0,
    "Y": 257.0,
    "Yb": 283.0,
    "Zn": 276.0,
    "Zr": 267.0
}

valence_dict = {
    "Ac": 11,
    "Ag": 11,
    "Al": 3,
    "Ar": 8,
    "As": 5,
    "At": 17,
    "Au": 11,
    "B": 3,
    "Ba": 10,
    "Be": 2,
    "Bi": 15,
    "Br": 7,
    "C": 4,
    "Ca": 10,
    "Cd": 12,
    "Ce": 11,
    "Cl": 7,
    "Co": 17,
    "Cr": 14,
    "Cs": 9,
    "Cu": 11,
    "Dy": 11,
    "Er": 11,
    "Eu": 11,
    "F": 7,
    "Fe": 16,
    "Fr": 19,
    "Ga": 13,
    "Gd": 11,
    "Ge": 14,
    "H": 1,
    "He": 2,
    "Hf": 12,
    "Hg": 12,
    "Ho": 11,
    "I": 7,
    "In": 13,
    "Ir": 9,
    "K": 9,
    "Kr": 18,
    "La": 11,
    "Li": 3,
    "Lu": 11,
    "Mg": 2,
    "Mn": 15,
    "Mo": 14,
    "N": 5,
    "Na": 9,
    "Nb": 13,
    "Nd": 11,
    "Ne": 8,
    "Ni": 10,
    "Np": 15,
    "O": 6,
    "Os": 16,
    "P": 5,
    "Pa": 13,
    "Pb": 14,
    "Pd": 10,
    "Pm": 11,
    "Po": 16,
    "Pr": 11,
    "Pt": 10,
    "Pu": 16,
    "Ra": 20,
    "Rb": 9,
    "Re": 15,
    "Rh": 17,
    "Rn": 18,
    "Ru": 16,
    "S": 6,
    "Sb": 5,
    "Sc": 11,
    "Se": 6,
    "Si": 4,
    "Sm": 11,
    "Sn": 14,
    "Sr": 10,
    "Ta": 13,
    "Tb": 11,
    "Tc": 15,
    "Te": 16,
    "Th": 12,
    "Ti": 12,
    "Tl": 13,
    "Tm": 11,
    "U": 14,
    "V": 13,
    "W": 14,
    "Xe": 18,
    "Y": 11,
    "Yb": 10,
    "Zn": 12,
    "Zr": 12
}

atomwfc_dict = {
    "Ac": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["7P", 3], ["6D", 5]]),
    "Ag": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Al": numpy.array([["3S", 1], ["3P", 3]]),
    "Ar": numpy.array([["3S", 1], ["3P", 3]]),
    "As": numpy.array([["4S", 1], ["4P", 3]]),
    "At": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Au": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "B": numpy.array([["2S", 1], ["2P", 3]]),
    "Ba": numpy.array([["5S", 1], ["6S", 1], ["5P", 3]]),
    "Be": numpy.array([["2S", 1], ["2P", 3]]),
    "Bi": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Br": numpy.array([["4S", 1], ["4P", 3]]),
    "C": numpy.array([["2S", 1], ["2P", 3]]),
    "Ca": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3]]),
    "Cd": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Ce": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Cl": numpy.array([["3S", 1], ["3P", 3]]),
    "Co": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "Cr": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "Cs": numpy.array([["5S", 1], ["6S", 1], ["5P", 3]]),
    "Cu": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Dy": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Er": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Eu": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "F": numpy.array([["2S", 1], ["2P", 3]]),
    "Fe": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3], ["3D", 5]]),
    "Fr": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["5D", 5]]),
    "Ga": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Gd": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Ge": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "H": numpy.array([["1S", 1]]),
    "He": numpy.array([["1S", 1]]),
    "Hf": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["5D", 5]]),
    "Hg": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Ho": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "I": numpy.array([["5S", 1], ["5P", 3]]),
    "In": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Ir": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "K": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3]]),
    "Kr": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "La": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["4F", 7]]),
    "Li": numpy.array([["1S", 1], ["2S", 1]]),
    "Lu": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Mg": numpy.array([["3S", 1], ["3P", 3]]),
    "Mn": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "Mo": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "N": numpy.array([["2S", 1], ["2P", 3]]),
    "Na": numpy.array([["2S", 1], ["3S", 1], ["2P", 3], ["3P", 3]]),
    "Nb": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Nd": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Ne": numpy.array([["2S", 1], ["2P", 3]]),
    "Ni": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Np": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "O": numpy.array([["2S", 1], ["2P", 3]]),
    "Os": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["5D", 5]]),
    "P": numpy.array([["3S", 1], ["3P", 3]]),
    "Pa": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "Pb": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Pd": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Pm": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Po": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Pr": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Pt": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Pu": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "Ra": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["5D", 5]]),
    "Rb": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3]]),
    "Re": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["5D", 5]]),
    "Rh": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Rn": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Ru": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["4D", 5]]),
    "S": numpy.array([["3S", 1], ["3P", 3]]),
    "Sb": numpy.array([["5S", 1], ["5P", 3]]),
    "Sc": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "Se": numpy.array([["4S", 1], ["4P", 3]]),
    "Si": numpy.array([["3S", 1], ["3P", 3]]),
    "Sm": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Sn": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Sr": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3]]),
    "Ta": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Tb": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Tc": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Te": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Th": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "Ti": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "Tl": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Tm": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "U": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "V": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "W": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Xe": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Y": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Yb": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Zn": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Zr": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]])
}

psen_dict = {
    "Ac": -9.37407640800e+01,
    "Ag": -8.40817940200e+01,
    "Al": -5.21990135000e+00,
    "Ar": -4.31314191200e+01,
    "As": -1.78981976300e+01,
    "At": -2.20357228920e+02,
    "Au": -8.51688880500e+01,
    "B": -5.81379974000e+00,
    "Ba": -6.25893335900e+01,
    "Be": -2.26997529000e+00,
    "Bi": -1.54936285780e+02,
    "Br": -3.37695477900e+01,
    "C": -1.13119748200e+01,
    "Ca": -7.49140824900e+01,
    "Cd": -1.16827249740e+02,
    "Ce": -1.26640727820e+02,
    "Cl": -3.45172024400e+01,
    "Co": -2.91099278750e+02,
    "Cr": -1.74392795640e+02,
    "Cu": -1.17148467640e+02,
    "Dy": -1.24634022130e+02,
    "Er": -1.22841427760e+02,
    "Eu": -1.23581941720e+02,
    "F": -4.92944940200e+01,
    "Fe": -2.53837149650e+02,
    "Fr": -2.84483981270e+02,
    "Ga": -1.72580703130e+02,
    "Gd": -1.23784180100e+02,
    "Ge": -2.11110250520e+02,
    "H": -9.31781210000e-01,
    "He": -5.75068890000e+00,
    "Hf": -1.24370465280e+02,
    "Hg": -1.06275525390e+02,
    "Ho": -1.22373632670e+02,
    "I": -5.88279169700e+01,
    "In": -1.44675602580e+02,
    "Ir": -6.14449998200e+01,
    "K": -5.87132036700e+01,
    "Kr": -3.42821144210e+02,
    "La": -1.25316377670e+02,
    "Li": -1.47244381400e+01,
    "Lu": -1.24708224560e+02,
    "Mg": -2.32745946000e+00,
    "Mn": -2.13271043220e+02,
    "Mo": -1.46979796590e+02,
    "N": -1.96190993500e+01,
    "Na": -9.21859679400e+01,
    "Nb": -1.17385335090e+02,
    "Nd": -1.25153999150e+02,
    "Ne": -6.99007359200e+01,
    "Ni": -1.00201571250e+02,
    "Np": -2.04716383260e+02,
    "O": -3.32066649000e+01,
    "Os": -1.98668898620e+02,
    "P": -1.51647862800e+01,
    "Pa": -1.14793811300e+02,
    "Pb": -1.46627668400e+02,
    "Pd": -7.84777686000e+01,
    "Pm": -1.25550497550e+02,
    "Po": -1.95388451910e+02,
    "Pr": -8.40940745800e+01,
    "Pt": -7.22671199700e+01,
    "Pu": -1.63557559670e+02,
    "Ra": -3.15536013120e+02,
    "Rb": -5.91270549700e+01,
    "Re": -1.77264265720e+02,
    "Rh": -2.32709435520e+02,
    "Rn": -2.37419043680e+02,
    "Ru": -1.96513393290e+02,
    "S": -2.25578742300e+01,
    "Sb": -1.96256771200e+01,
    "Sc": -9.75730792100e+01,
    "Se": -2.74714561300e+01,
    "Si": -1.10573775700e+01,
    "Sm": -1.25928165630e+02,
    "Sn": -1.48720541990e+02,
    "Sr": -7.00672433000e+01,
    "Ta": -1.52548382960e+02,
    "Tb": -1.24147531950e+02,
    "Tc": -1.79570929440e+02,
    "Te": -2.14917232170e+02,
    "Th": -1.55389683780e+02,
    "Ti": -1.19293928420e+02,
    "Tl": -1.16696457060e+02,
    "Tm": -1.23487328960e+02,
    "U": -1.28820990520e+02,
    "V": -1.45000178730e+02,
    "W": -1.57147131880e+02,
    "Xe": -2.82239895140e+02,
    "Y": -8.31798527300e+01,
    "Yb": -1.13761514660e+02,
    "Zn": -1.48256520750e+02,
    "Zr": -9.89062225400e+01
}

band_dict = {
    "Ac": 11,
    "Ag": 6,
    "Al": 4,
    "Ar": 4,
    "As": 4,
    "At": 17,
    "Au": 6,
    "B": 4,
    "Ba": 8,
    "Be": 4,
    "Bi": 15,
    "Br": 4,
    "C": 4,
    "Ca": 8,
    "Cd": 12,
    "Ce": 11,
    "Cl": 4,
    "Co": 17,
    "Cr": 14,
    "Cs": 8,
    "Cu": 11,
    "Dy": 11,
    "Er": 11,
    "Eu": 11,
    "F": 4,
    "Fe": 16,
    "Fr": 19,
    "Ga": 13,
    "Gd": 11,
    "Ge": 14,
    "H": 4,
    "He": 1,
    "Hf": 12,
    "Hg": 12,
    "Ho": 11,
    "I": 4,
    "In": 13,
    "Ir": 9,
    "K": 8,
    "Kr": 9,
    "La": 11,
    "Li": 5,
    "Lu": 11,
    "Mg": 4,
    "Mn": 10,
    "Mo": 10,
    "N": 4,
    "Na": 9,
    "Nb": 10,
    "Nd": 11,
    "Ne": 8,
    "Ni": 6,
    "Np": 15,
    "O": 4,
    "Os": 16,
    "P": 4,
    "Pa": 13,
    "Pb": 14,
    "Pd": 6,
    "Pm": 11,
    "Po": 16,
    "Pr": 11,
    "Pt": 10,
    "Pu": 16,
    "Ra": 20,
    "Rb": 9,
    "Re": 15,
    "Rh": 17,
    "Rn": 18,
    "Ru": 16,
    "S": 4,
    "Sb": 5,
    "Sc": 11,
    "Se": 4,
    "Si": 4,
    "Sm": 11,
    "Sn": 14,
    "Sr": 10,
    "Ta": 13,
    "Tb": 11,
    "Tc": 15,
    "Te": 16,
    "Th": 12,
    "Ti": 12,
    "Tl": 13,
    "Tm": 11,
    "U": 14,
    "V": 13,
    "W": 14,
    "Xe": 18,
    "Y": 11,
    "Yb": 10,
    "Zn": 12,
    "Zr": 12
}
