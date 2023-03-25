import numpy


pseudo_dict = {
    "B": "B.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "C": "C.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "F": "F.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "H": "H.rel-pbe-rrkjus_psl.0.1.UPF",
    "I": "I.rel-pbe-n-rrkjus_psl.0.2.2.UPF",
    "K": "K.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "N": "N.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "O": "O.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "P": "P.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "S": "S.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "U": "U.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "V": "V.rel-pbe-spnl-rrkjus_psl.1.0.0.UPF",
    "W": "W.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Y": "Y.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Ac": "Ac.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Ag": "Ag.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
    "Al": "Al.rel-pbe-n-rrkjus_psl.0.2.2.UPF",
    "Ar": "Ar.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
    "As": "As.rel-pbe-n-rrkjus_psl.0.2.UPF",
    "At": "At.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Au": "Au.rel-pbe-dn-rrkjus_psl.0.1.UPF",
    "Ba": "Ba.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Be": "Be.rel-pbe-n-rrkjus_psl.0.2.UPF",
    "Bi": "Bi.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Br": "Br.rel-pbe-n-rrkjus_psl.0.2.UPF",
    "Ca": "Ca.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Cd": "Cd.rel-pbe-dn-rrkjus_psl.0.2.UPF",
    "Ce": "Ce.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Cl": "Cl.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "Co": "Co.rel-pbe-spn-rrkjus_psl.0.3.1.UPF",
    "Cr": "Cr.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Cs": "Cs.rel-pbe-spnl-rrkjus_psl.1.0.0.UPF",
    "Cu": "Cu.rel-pbe-dn-rrkjus_psl.0.2.UPF",
    "Dy": "Dy.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Er": "Er.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Eu": "Eu.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Fe": "Fe.rel-pbe-spn-rrkjus_psl.0.2.1.UPF",
    "Fr": "Fr.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Ga": "Ga.rel-pbe-dn-rrkjus_psl.0.2.UPF",
    "Gd": "Gd.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Ge": "Ge.rel-pbe-dn-rrkjus_psl.0.2.2.UPF",
    "He": "He.rel-pbe-rrkjus_psl.1.0.0.UPF",
    "Hf": "Hf.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Hg": "Hg.rel-pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Ho": "Ho.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "In": "In.rel-pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Ir": "Ir.rel-pbe-n-rrkjus_psl.0.2.3.UPF",
    "Kr": "Kr.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "La": "La.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Li": "Li.rel-pbe-s-rrkjus_psl.0.2.1.UPF",
    "Lu": "Lu.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Mg": "Mg.rel-pbe-n-rrkjus_psl.0.3.0.UPF",
    "Mn": "Mn.rel-pbe-spn-rrkjus_psl.0.3.1.UPF",
    "Mo": "Mo.rel-pbe-spn-rrkjus_psl.0.2.UPF",
    "Na": "Na.rel-pbe-spn-rrkjus_psl.0.2.UPF",
    "Nb": "Nb.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Nd": "Nd.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Ne": "Ne.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
    "Ni": "Ni.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "Np": "Np.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Os": "Os.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Pa": "Pa.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Pb": "Pb.rel-pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Pd": "Pd.rel-pbe-n-rrkjus_psl.0.2.2.UPF",
    "Pm": "Pm.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Po": "Po.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Pr": "Pr.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Pt": "Pt.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "Pu": "Pu.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Ra": "Ra.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Rb": "Rb.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Re": "Re.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Rh": "Rh.rel-pbe-spn-rrkjus_psl.0.2.3.UPF",
    "Rn": "Rn.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Ru": "Ru.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Sb": "Sb.rel-pbe-n-rrkjus_psl.1.0.0.UPF",
    "Sc": "Sc.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Se": "Se.rel-pbe-n-rrkjus_psl.0.2.UPF",
    "Si": "Si.rel-pbe-n-rrkjus_psl.0.1.UPF",
    "Sm": "Sm.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Sn": "Sn.rel-pbe-dn-rrkjus_psl.0.2.UPF",
    "Sr": "Sr.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Ta": "Ta.rel-pbe-spn-rrkjus_psl.0.2.UPF",
    "Tb": "Tb.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Tc": "Tc.rel-pbe-spn-rrkjus_psl.0.3.0.UPF",
    "Te": "Te.rel-pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Th": "Th.rel-pbe-spfn-rrkjus_psl.1.0.0.UPF",
    "Ti": "Ti.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Tl": "Tl.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Tm": "Tm.rel-pbe-spdn-rrkjus_psl.1.0.0.UPF",
    "Xe": "Xe.rel-pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Yb": "Yb.rel-pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Zn": "Zn.rel-pbe-dnl-rrkjus_psl.1.0.0.UPF",
    "Zr": "Zr.rel-pbe-spn-rrkjus_psl.1.0.0.UPF"
    }

ecutwfc_dict = {
    "B": 48.0,
    "C": 37.0,
    "F": 48.0,
    "H": 46.0,
    "I": 28.0,
    "K": 41.0,
    "N": 39.0,
    "O": 47.0,
    "P": 17.0,
    "S": 18.0,
    "U": 93.0,
    "V": 48.0,
    "W": 55.0,
    "Y": 39.0,
    "Ac": 48.0,
    "Ag": 45.0,
    "Al": 37.0,
    "Ar": 48.0,
    "As": 20.0,
    "At": 50.0,
    "Au": 27.0,
    "Ba": 38.0,
    "Be": 30.0,
    "Bi": 45.0,
    "Br": 25.0,
    "Ca": 45.0,
    "Cd": 31.0,
    "Ce": 37.0,
    "Cl": 36.0,
    "Co": 60.0,
    "Cr": 48.0,
    "Cs": 35.0,
    "Cu": 33.0,
    "Dy": 38.0,
    "Er": 41.0,
    "Eu": 35.0,
    "Fe": 64.0,
    "Fr": 67.0,
    "Ga": 33.0,
    "Gd": 35.0,
    "Ge": 38.0,
    "He": 46.0,
    "Hf": 47.0,
    "Hg": 30.0,
    "Ho": 39.0,
    "In": 48.0,
    "Ir": 40.0,
    "Kr": 55.0,
    "La": 75.0,
    "Li": 78.0,
    "Lu": 46.0,
    "Mg": 13.0,
    "Mn": 46.0,
    "Mo": 49.0,
    "Na": 41.0,
    "Nb": 48.0,
    "Nd": 38.0,
    "Ne": 54.0,
    "Ni": 41.0,
    "Np": 111.0,
    "Os": 54.0,
    "Pa": 91.0,
    "Pb": 41.0,
    "Pd": 32.0,
    "Pm": 34.0,
    "Po": 85.0,
    "Pr": 41.0,
    "Pt": 40.0,
    "Pu": 45.0,
    "Ra": 82.0,
    "Rb": 35.0,
    "Re": 54.0,
    "Rh": 40.0,
    "Rn": 50.0,
    "Ru": 52.0,
    "Sb": 27.0,
    "Sc": 50.0,
    "Se": 21.0,
    "Si": 38.0,
    "Sm": 34.0,
    "Sn": 25.0,
    "Sr": 39.0,
    "Ta": 53.0,
    "Tb": 37.0,
    "Tc": 63.0,
    "Te": 34.0,
    "Th": 103.0,
    "Ti": 52.0,
    "Tl": 45.0,
    "Tm": 42.0,
    "Xe": 45.0,
    "Yb": 41.0,
    "Zn": 44.0,
    "Zr": 46.0
}

ecutrho_dict = {
    "B": 216.0,
    "C": 318.0,
    "F": 316.0,
    "H": 221.0,
    "I": 113.0,
    "K": 277.0,
    "N": 264.0,
    "O": 323.0,
    "P": 79.0,
    "S": 77.0,
    "U": 701.0,
    "V": 645.0,
    "W": 373.0,
    "Y": 257.0,
    "Ac": 329.0,
    "Ag": 182.0,
    "Al": 146.0,
    "Ar": 225.0,
    "As": 103.0,
    "At": 473.0,
    "Au": 391.0,
    "Ba": 192.0,
    "Be": 197.0,
    "Bi": 456.0,
    "Br": 108.0,
    "Ca": 274.0,
    "Cd": 156.0,
    "Ce": 179.0,
    "Cl": 146.0,
    "Co": 450.0,
    "Cr": 457.0,
    "Cs": 173.0,
    "Cu": 150.0,
    "Dy": 266.0,
    "Er": 267.0,
    "Eu": 259.0,
    "Fe": 782.0,
    "Fr": 508.0,
    "Ga": 180.0,
    "Gd": 227.0,
    "Ge": 238.0,
    "He": 203.0,
    "Hf": 265.0,
    "Hg": 121.0,
    "Ho": 269.0,
    "In": 193.0,
    "Ir": 526.0,
    "Kr": 252.0,
    "La": 537.0,
    "Li": 355.0,
    "Lu": 346.0,
    "Mg": 87.0,
    "Mn": 244.0,
    "Mo": 408.0,
    "Na": 163.0,
    "Nb": 266.0,
    "Nd": 202.0,
    "Ne": 264.0,
    "Ni": 236.0,
    "Np": 561.0,
    "Os": 270.0,
    "Pa": 447.0,
    "Pb": 166.0,
    "Pd": 212.0,
    "Pm": 198.0,
    "Po": 457.0,
    "Pr": 310.0,
    "Pt": 433.0,
    "Pu": 353.0,
    "Ra": 534.0,
    "Rb": 259.0,
    "Re": 290.0,
    "Rh": 439.0,
    "Rn": 220.0,
    "Ru": 368.0,
    "Sb": 242.0,
    "Sc": 402.0,
    "Se": 95.0,
    "Si": 151.0,
    "Sm": 195.0,
    "Sn": 137.0,
    "Sr": 263.0,
    "Ta": 402.0,
    "Tb": 229.0,
    "Tc": 853.0,
    "Te": 148.0,
    "Th": 414.0,
    "Ti": 575.0,
    "Tl": 210.0,
    "Tm": 319.0,
    "Xe": 212.0,
    "Yb": 283.0,
    "Zn": 276.0,
    "Zr": 268.0
}

valence_dict = {
    "B": 3,
    "C": 4,
    "F": 7,
    "H": 1,
    "I": 7,
    "K": 9,
    "N": 5,
    "O": 6,
    "P": 5,
    "S": 6,
    "U": 14,
    "V": 13,
    "W": 14,
    "Y": 11,
    "Ac": 11,
    "Ag": 11,
    "Al": 3,
    "Ar": 8,
    "As": 5,
    "At": 17,
    "Au": 11,
    "Ba": 10,
    "Be": 2,
    "Bi": 15,
    "Br": 7,
    "Ca": 10,
    "Cd": 12,
    "Ce": 11,
    "Cl": 7,
    "Co": 17,
    "Cr": 14,
    "Cs": -5,
    "Cu": 11,
    "Dy": 11,
    "Er": 11,
    "Eu": 11,
    "Fe": 16,
    "Fr": 19,
    "Ga": 13,
    "Gd": 11,
    "Ge": 14,
    "He": 2,
    "Hf": 12,
    "Hg": 12,
    "Ho": 11,
    "In": 13,
    "Ir": 9,
    "Kr": 18,
    "La": 11,
    "Li": 3,
    "Lu": 11,
    "Mg": 2,
    "Mn": 15,
    "Mo": 14,
    "Na": 9,
    "Nb": 13,
    "Nd": 11,
    "Ne": 8,
    "Ni": 10,
    "Np": 15,
    "Os": 16,
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
    "Xe": 18,
    "Yb": 10,
    "Zn": 12,
    "Zr": 12
}

atomwfc_dict = {
    "B": numpy.array([["2S", 1], ["2P", 3]]),
    "C": numpy.array([["2S", 1], ["2P", 3]]),
    "F": numpy.array([["2S", 1], ["2P", 3]]),
    "H": numpy.array([["1S", 1]]),
    "I": numpy.array([["5S", 1], ["5P", 3]]),
    "K": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3]]),
    "N": numpy.array([["2S", 1], ["2P", 3]]),
    "O": numpy.array([["2S", 1], ["2P", 3]]),
    "P": numpy.array([["3S", 1], ["3P", 3]]),
    "S": numpy.array([["3S", 1], ["3P", 3]]),
    "U": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "V": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "W": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Y": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Ac": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["7P", 3], ["6D", 5]]),
    "Ag": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Al": numpy.array([["3S", 1], ["3P", 3]]),
    "Ar": numpy.array([["3S", 1], ["3P", 3]]),
    "As": numpy.array([["4S", 1], ["4P", 3]]),
    "At": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Au": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Ba": numpy.array([["5S", 1], ["6S", 1], ["5P", 3]]),
    "Be": numpy.array([["2S", 1], ["2P", 3]]),
    "Bi": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Br": numpy.array([["4S", 1], ["4P", 3]]),
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
    "Fe": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3], ["3D", 5]]),
    "Fr": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["5D", 5]]),
    "Ga": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Gd": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Ge": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "He": numpy.array([["1S", 1]]),
    "Hf": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["5D", 5]]),
    "Hg": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Ho": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "In": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Ir": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Kr": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "La": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["4F", 7]]),
    "Li": numpy.array([["1S", 1], ["2S", 1]]),
    "Lu": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Mg": numpy.array([["3S", 1], ["3P", 3]]),
    "Mn": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["3D", 5]]),
    "Mo": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Na": numpy.array([["2S", 1], ["3S", 1], ["2P", 3], ["3P", 3]]),
    "Nb": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]]),
    "Nd": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Ne": numpy.array([["2S", 1], ["2P", 3]]),
    "Ni": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Np": numpy.array([["6S", 1], ["7S", 1], ["6P", 3], ["6D", 5], ["5F", 7]]),
    "Os": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["5D", 5]]),
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
    "Xe": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Yb": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5]]),
    "Zn": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "Zr": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["5P", 3], ["4D", 5]])
}
