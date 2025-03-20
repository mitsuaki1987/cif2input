import numpy


pseudo_dict = {
    "Ag": "Ag_ONCV_PBE-1.0.oncvpsp.upf",
    "Al": "Al.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Ar": "Ar_ONCV_PBE-1.1.oncvpsp.upf",
    "As": "As.pbe-n-rrkjus_psl.0.2.UPF",
    "Au": "Au_ONCV_PBE-1.0.oncvpsp.upf",
    "B": "b_pbe_v1.4.uspp.F.UPF",
    "Ba": "Ba.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Be": "be_pbe_v1.4.uspp.F.UPF",
    "Bi": "Bi_pbe_v1.uspp.F.UPF",
    "Br": "br_pbe_v1.4.uspp.F.UPF",
    "C": "C.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Ca": "Ca_pbe_v1.uspp.F.UPF",
    "Cd": "Cd.pbe-dn-rrkjus_psl.0.3.1.UPF",
    "Ce": "Ce.GGA-PBE-paw-v1.0.UPF",
    "Cl": "cl_pbe_v1.4.uspp.F.UPF",
    "Co": "Co_pbe_v1.2.uspp.F.UPF",
    "Cr": "cr_pbe_v1.5.uspp.F.UPF",
    "Cs": "Cs_pbe_v1.uspp.F.UPF",
    "Cu": "Cu_pbe_v1.2.uspp.F.UPF",
    "Dy": "Dy.GGA-PBE-paw-v1.0.UPF",
    "Er": "Er.GGA-PBE-paw-v1.0.UPF",
    "Eu": "Eu.GGA-PBE-paw-v1.0.UPF",
    "F": "f_pbe_v1.4.uspp.F.UPF",
    "Fe": "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF",
    "Ga": "Ga.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Gd": "Gd.GGA-PBE-paw-v1.0.UPF",
    "Ge": "ge_pbe_v1.4.uspp.F.UPF",
    "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
    "He": "He_ONCV_PBE-1.0.oncvpsp.upf",
    "Hf": "Hf-sp.oncvpsp.upf",
    "Hg": "Hg_ONCV_PBE-1.0.oncvpsp.upf",
    "Ho": "Ho.GGA-PBE-paw-v1.0.UPF",
    "I": "I.pbe-n-rrkjus_psl.0.2.UPF",
    "In": "In.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Ir": "Ir_pbe_v1.2.uspp.F.UPF",
    "K": "K.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Kr": "Kr_ONCV_PBE-1.0.oncvpsp.upf",
    "La": "La.GGA-PBE-paw-v1.0.UPF",
    "Li": "li_pbe_v1.4.uspp.F.UPF",
    "Lu": "Lu.GGA-PBE-paw-v1.0.UPF",
    "Mg": "Mg.pbe-n-rrkjus_psl.0.3.0.UPF",
    "Mn": "mn_pbe_v1.5.uspp.F.UPF",
    "Mo": "Mo_ONCV_PBE-1.0.oncvpsp.upf",
    "N": "N.pbe-n-radius_5.UPF",
    "Na": "na_pbe_v1.5.uspp.F.UPF",
    "Nb": "Nb.pbe-spn-rrkjus_psl.0.3.0.UPF",
    "Nd": "Nd.GGA-PBE-paw-v1.0.UPF",
    "Ne": "Ne_ONCV_PBE-1.0.oncvpsp.upf",
    "Ni": "ni_pbe_v1.4.uspp.F.UPF",
    "O": "O.pbe-n-rrkjus_psl.0.1.UPF",
    "Os": "Os_pbe_v1.2.uspp.F.UPF",
    "P": "P.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Pb": "Pb.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Pd": "Pd_ONCV_PBE-1.0.oncvpsp.upf",
    "Pm": "Pm.GGA-PBE-paw-v1.0.UPF",
    "Po": "Po.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Pr": "Pr.GGA-PBE-paw-v1.0.UPF",
    "Pt": "pt_pbe_v1.4.uspp.F.UPF",
    "Rb": "Rb_ONCV_PBE-1.0.oncvpsp.upf",
    "Re": "Re_pbe_v1.2.uspp.F.UPF",
    "Rh": "Rh_ONCV_PBE-1.0.oncvpsp.upf",
    "Rn": "Rn.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "Ru": "Ru_ONCV_PBE-1.0.oncvpsp.upf",
    "S": "s_pbe_v1.4.uspp.F.UPF",
    "Sb": "sb_pbe_v1.4.uspp.F.UPF",
    "Sc": "Sc_ONCV_PBE-1.0.oncvpsp.upf",
    "Se": "Se_pbe_v1.uspp.F.UPF",
    "Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Sm": "Sm.GGA-PBE-paw-v1.0.UPF",
    "Sn": "Sn_pbe_v1.uspp.F.UPF",
    "Sr": "Sr_pbe_v1.uspp.F.UPF",
    "Ta": "Ta_pbe_v1.uspp.F.UPF",
    "Tb": "Tb.GGA-PBE-paw-v1.0.UPF",
    "Tc": "Tc_ONCV_PBE-1.0.oncvpsp.upf",
    "Te": "Te_pbe_v1.uspp.F.UPF",
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "Tl": "Tl_pbe_v1.2.uspp.F.UPF",
    "Tm": "Tm.GGA-PBE-paw-v1.0.UPF",
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "W": "W_pbe_v1.2.uspp.F.UPF",
    "Xe": "Xe_ONCV_PBE-1.1.oncvpsp.upf",
    "Y": "Y_pbe_v1.uspp.F.UPF",
    "Yb": "Yb.GGA-PBE-paw-v1.0.UPF",
    "Zn": "Zn_pbe_v1.uspp.F.UPF",
    "Zr": "Zr_pbe_v1.uspp.F.UPF"
}

ecutwfc_dict = {
    "Ag": 50.0,
    "Al": 30.0,
    "Ar": 60.0,
    "As": 35.0,
    "Au": 45.0,
    "B": 35.0,
    "Ba": 30.0,
    "Be": 40.0,
    "Bi": 45.0,
    "Br": 30.0,
    "C": 45.0,
    "Ca": 30.0,
    "Cd": 60.0,
    "Ce": 40.0,
    "Cl": 40.0,
    "Co": 45.0,
    "Cr": 40.0,
    "Cs": 30.0,
    "Cu": 55.0,
    "Dy": 40.0,
    "Er": 40.0,
    "Eu": 40.0,
    "F": 45.0,
    "Fe": 90.0,
    "Ga": 70.0,
    "Gd": 40.0,
    "Ge": 40.0,
    "H": 60.0,
    "He": 50.0,
    "Hf": 50.0,
    "Hg": 50.0,
    "Ho": 40.0,
    "I": 35.0,
    "In": 50.0,
    "Ir": 55.0,
    "K": 60.0,
    "Kr": 45.0,
    "La": 40.0,
    "Li": 40.0,
    "Lu": 45.0,
    "Mg": 30.0,
    "Mn": 65.0,
    "Mo": 35.0,
    "N": 60.0,
    "Na": 40.0,
    "Nb": 40.0,
    "Nd": 40.0,
    "Ne": 50.0,
    "Ni": 45.0,
    "O": 50.0,
    "Os": 40.0,
    "P": 30.0,
    "Pb": 40.0,
    "Pd": 45.0,
    "Pm": 40.0,
    "Po": 75.0,
    "Pr": 40.0,
    "Pt": 35.0,
    "Rb": 30.0,
    "Re": 30.0,
    "Rh": 35.0,
    "Rn": 120.0,
    "Ru": 35.0,
    "S": 35.0,
    "Sb": 40.0,
    "Sc": 40.0,
    "Se": 30.0,
    "Si": 30.0,
    "Sm": 40.0,
    "Sn": 60.0,
    "Sr": 30.0,
    "Ta": 45.0,
    "Tb": 40.0,
    "Tc": 30.0,
    "Te": 30.0,
    "Ti": 35.0,
    "Tl": 50.0,
    "Tm": 40.0,
    "V": 35.0,
    "W": 30.0,
    "Xe": 60.0,
    "Y": 35.0,
    "Yb": 40.0,
    "Zn": 40.0,
    "Zr": 30.0
}

ecutrho_dict = {
    "Ag": 200.0,
    "Al": 240.0,
    "Ar": 240.0,
    "As": 280.0,
    "Au": 180.0,
    "B": 280.0,
    "Ba": 240.0,
    "Be": 320.0,
    "Bi": 360.0,
    "Br": 240.0,
    "C": 360.0,
    "Ca": 240.0,
    "Cd": 480.0,
    "Ce": 320.0,
    "Cl": 320.0,
    "Co": 360.0,
    "Cr": 320.0,
    "Cs": 240.0,
    "Cu": 440.0,
    "Dy": 320.0,
    "Er": 320.0,
    "Eu": 320.0,
    "F": 360.0,
    "Fe": 1080.0,
    "Ga": 560.0,
    "Gd": 320.0,
    "Ge": 320.0,
    "H": 480.0,
    "He": 200.0,
    "Hf": 200.0,
    "Hg": 200.0,
    "Ho": 320.0,
    "I": 280.0,
    "In": 400.0,
    "Ir": 440.0,
    "K": 480.0,
    "Kr": 180.0,
    "La": 320.0,
    "Li": 320.0,
    "Lu": 360.0,
    "Mg": 240.0,
    "Mn": 780.0,
    "Mo": 140.0,
    "N": 480.0,
    "Na": 320.0,
    "Nb": 320.0,
    "Nd": 320.0,
    "Ne": 200.0,
    "Ni": 360.0,
    "O": 800.0,
    "Os": 320.0,
    "P": 240.0,
    "Pb": 320.0,
    "Pd": 180.0,
    "Pm": 320.0,
    "Po": 600.0,
    "Pr": 320.0,
    "Pt": 280.0,
    "Rb": 120.0,
    "Re": 240.0,
    "Rh": 140.0,
    "Rn": 960.0,
    "Ru": 140.0,
    "S": 280.0,
    "Sb": 320.0,
    "Sc": 160.0,
    "Se": 240.0,
    "Si": 240.0,
    "Sm": 320.0,
    "Sn": 480.0,
    "Sr": 240.0,
    "Ta": 360.0,
    "Tb": 320.0,
    "Tc": 120.0,
    "Te": 240.0,
    "Ti": 280.0,
    "Tl": 400.0,
    "Tm": 320.0,
    "V": 280.0,
    "W": 240.0,
    "Xe": 240.0,
    "Y": 280.0,
    "Yb": 320.0,
    "Zn": 320.0,
    "Zr": 240.0
}

valence_dict = {
    "Ag": 19,
    "Al": 3,
    "Ar": 8,
    "As": 5,
    "Au": 19,
    "B": 3,
    "Ba": 10,
    "Be": 4,
    "Bi": 15,
    "Br": 7,
    "C": 4,
    "Ca": 10,
    "Cd": 12,
    "Ce": 12,
    "Cl": 7,
    "Co": 17,
    "Cr": 14,
    "Cs": 9,
    "Cu": 19,
    "Dy": 20,
    "Er": 22,
    "Eu": 17,
    "F": 7,
    "Fe": 16,
    "Ga": 13,
    "Gd": 18,
    "Ge": 14,
    "H": 1,
    "He": 2,
    "Hf": 12,
    "Hg": 20,
    "Ho": 21,
    "I": 7,
    "In": 13,
    "Ir": 15,
    "K": 9,
    "Kr": 8,
    "La": 11,
    "Li": 3,
    "Lu": 25,
    "Mg": 2,
    "Mn": 15,
    "Mo": 14,
    "N": 5,
    "Na": 9,
    "Nb": 13,
    "Nd": 14,
    "Ne": 8,
    "Ni": 18,
    "O": 6,
    "Os": 16,
    "P": 5,
    "Pb": 14,
    "Pd": 18,
    "Pm": 15,
    "Po": 16,
    "Pr": 13,
    "Pt": 16,
    "Rb": 9,
    "Re": 15,
    "Rh": 17,
    "Rn": 18,
    "Ru": 16,
    "S": 6,
    "Sb": 15,
    "Sc": 11,
    "Se": 6,
    "Si": 4,
    "Sm": 16,
    "Sn": 14,
    "Sr": 10,
    "Ta": 13,
    "Tb": 19,
    "Tc": 15,
    "Te": 6,
    "Ti": 12,
    "Tl": 13,
    "Tm": 23,
    "V": 13,
    "W": 14,
    "Xe": 18,
    "Y": 11,
    "Yb": 24,
    "Zn": 20,
    "Zr": 12
}

atomwfc_dict = {
    "N": numpy.array([["2S", 1], ["2P", 3]]),
    "Ce": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Dy": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Er": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Eu": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Gd": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Ho": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "La": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Lu": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Nd": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Pm": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Pr": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Sm": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Tb": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Tm": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Yb": numpy.array([["5S", 1], ["6S", 1], ["5P", 3], ["6P", 3], ["5D", 5], ["6D", 5], ["4F", 7], ["5F", 7]]),
    "Bi": numpy.array([["5D", 5], ["6S", 1], ["6P", 3]]),
    "Ca": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Co": numpy.array([["3S", 1], ["3P", 3], ["3D", 5], ["4S", 1], ["4P", 3]]),
    "Cs": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "Cu": numpy.array([["3S", 1], ["3P", 3], ["3D", 5], ["4S", 1], ["4P", 3]]),
    "Ir": numpy.array([["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "Os": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "Re": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "Se": numpy.array([["4S", 1], ["4P", 3]]),
    "Sn": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Sr": numpy.array([["4S", 1], ["4P", 3], ["4D", 5], ["5S", 1], ["5P", 3]]),
    "Ta": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "Te": numpy.array([["5S", 1], ["5P", 3]]),
    "Tl": numpy.array([["5D", 5], ["6S", 1], ["6P", 3]]),
    "W": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "Y": numpy.array([["4S", 1], ["4P", 3], ["4D", 5], ["5S", 1], ["5P", 3]]),
    "Zn": numpy.array([["3S", 1], ["3P", 3], ["3D", 5], ["4S", 1], ["4P", 3]]),
    "Zr": numpy.array([["4S", 1], ["4P", 3], ["4D", 5], ["5S", 1], ["5P", 3]]),
    "B": numpy.array([["2S", 1], ["2P", 3]]),
    "Be": numpy.array([["1S", 1], ["2S", 1], ["2P", 3]]),
    "Br": numpy.array([["4S", 1], ["4P", 3]]),
    "Cl": numpy.array([["3S", 1], ["3P", 3]]),
    "Cr": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "F": numpy.array([["2S", 1], ["2P", 3]]),
    "Ge": numpy.array([["3D", 5], ["4S", 1], ["4P", 3]]),
    "Li": numpy.array([["1S", 1], ["2S", 1], ["2P", 3]]),
    "Mn": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5], ["4P", 3]]),
    "Na": numpy.array([["2S", 1], ["2P", 3], ["3S", 1]]),
    "Ni": numpy.array([["3S", 1], ["3P", 3], ["3D", 5], ["4S", 1], ["4P", 3]]),
    "Pt": numpy.array([["5P", 3], ["5D", 5], ["6S", 1], ["6P", 3]]),
    "S": numpy.array([["3S", 1], ["3P", 3]]),
    "Sb": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Ti": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "V": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Hf": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1]]),
    "Ag": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Ar": numpy.array([["3S", 1], ["3P", 3]]),
    "Au": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "He": numpy.array([["1S", 1]]),
    "Hg": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "Kr": numpy.array([["4S", 1], ["4P", 3]]),
    "Mo": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Ne": numpy.array([["2S", 1], ["2P", 3]]),
    "Pd": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Rb": numpy.array([["4S", 1], ["4P", 3], ["5S", 1]]),
    "Rh": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Ru": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Sc": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Tc": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Xe": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Al": numpy.array([["3S", 1], ["3P", 3]]),
    "As": numpy.array([["4S", 1], ["4P", 3]]),
    "Ba": numpy.array([["5S", 1], ["6S", 1], ["5P", 3]]),
    "C": numpy.array([["2S", 1], ["2P", 3]]),
    "Cd": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "Fe": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3], ["3D", 5]]),
    "Ga": numpy.array([["4S", 1], ["4P", 3], ["3D", 5]]),
    "H": numpy.array([["1S", 1]]),
    "I": numpy.array([["5S", 1], ["5P", 3]]),
    "In": numpy.array([["5S", 1], ["5P", 3], ["4D", 5]]),
    "K": numpy.array([["3S", 1], ["4S", 1], ["3P", 3], ["4P", 3]]),
    "Mg": numpy.array([["3S", 1], ["3P", 3]]),
    "Nb": numpy.array([["4S", 1], ["5S", 1], ["4P", 3], ["4D", 5]]),
    "O": numpy.array([["2S", 1], ["2P", 3]]),
    "P": numpy.array([["3S", 1], ["3P", 3]]),
    "Pb": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Po": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Rn": numpy.array([["6S", 1], ["6P", 3], ["5D", 5]]),
    "Si": numpy.array([["3S", 1], ["3P", 3]]),
}

psen_dict = {
    "Ag": -2.87458497570e+02,
    "Al": -3.92485924200e+01,
    "Ar": -4.23888901800e+01,
    "As": -1.78986981500e+01,
    "Au": -2.74500089190e+02,
    "B": -5.92124151000e+00,
    "Ba": -5.13112242500e+01,
    "Be": -2.79694165200e+01,
    "Bi": -1.85030674010e+02,
    "Br": -4.67488004500e+01,
    "C": -1.78024970700e+01,
    "Ca": -7.47727607100e+01,
    "Cd": -1.22183071940e+02,
    "Ce": -4.81677998310e+02,
    "Cl": -3.31797287900e+01,
    "Co": -2.98022482630e+02,
    "Cr": -1.74895527070e+02,
    "Cs": -6.30187768000e+01,
    "Cu": -4.03508309890e+02,
    "Dy": -7.76876713810e+02,
    "Er": -8.85685358530e+02,
    "F": -4.83467469900e+01,
    "Fe": -3.28739479760e+02,
    "Ga": -2.77734401910e+02,
    "Gd": -6.83319769990e+02,
    "Ge": -2.13355823810e+02,
    "H": -9.31863260000e-01,
    "He": -5.58120348000e+00,
    "Hf": -1.12333451900e+02,
    "Hg": -3.10268913870e+02,
    "Ho": -8.29298947440e+02,
    "I": -3.78823193720e+02,
    "In": -1.44675700230e+02,
    "Ir": -1.81189080220e+02,
    "K": -1.12951344940e+02,
    "Kr": -3.71776222000e+01,
    "La": -4.57885292570e+02,
    "Li": -1.43471665500e+01,
    "Lu": -1.07985414802e+03,
    "Mg": -3.34977475200e+01,
    "Mn": -2.10636289760e+02,
    "Mo": -1.36441981150e+02,
    "N": -1.96183421200e+01,
    "Na": -9.52648947100e+01,
    "Nb": -3.33620624290e+02,
    "Nd": -5.36994239450e+02,
    "Ne": -6.13601036800e+01,
    "Ni": -3.42908667660e+02,
    "O": -4.11929764100e+01,
    "Os": -8.85018911080e+02,
    "P": -1.36933580400e+01,
    "Pb": -1.12671379520e+02,
    "Pd": -2.53104481200e+02,
    "Pm": -5.68866448330e+02,
    "Po": -1.95399311700e+02,
    "Pr": -5.07987931260e+02,
    "Pt": -2.10277583880e+02,
    "Rb": -4.88737685900e+01,
    "Re": -1.89873667920e+02,
    "Rh": -2.16669801870e+02,
    "Rn": -1.02708446964e+03,
    "Ru": -1.88069238110e+02,
    "S": -2.38341781100e+01,
    "Sb": -1.43593488140e+02,
    "Sc": -9.15689769100e+01,
    "Se": -4.33824584500e+01,
    "Si": -1.10591471900e+01,
    "Sm": -6.03784071880e+02,
    "Sn": -1.62997309110e+02,
    "Sr": -7.00151316400e+01,
    "Ta": -1.41223234700e+02,
    "Tb": -7.28269334970e+02,
    "Tc": -1.60575133670e+02,
    "Te": -2.63057042900e+01,
    "Ti": -1.18801564050e+02,
    "Tl": -1.44823534890e+02,
    "V": -1.44527758160e+02,
    "W": -1.58354378450e+02,
    "Xe": -2.31787053180e+02,
    "Y": -8.71802309400e+01,
    "Yb": -1.01093451642e+03,
    "Zn": -4.61462677840e+02,
    "Zr": -9.86560719800e+01
}

band_dict = {
    "Ag": 10,
    "Al": 4,
    "Ar": 4,
    "As": 4,
    "Au": 10,
    "B": 4,
    "Ba": 8,
    "Be": 5,
    "Bi": 9,
    "Br": 4,
    "C": 4,
    "Ca": 8,
    "Cd": 9,
    "Ce": 17,
    "Cl": 4,
    "Co": 10,
    "Cr": 10,
    "Cs": 8,
    "Cu": 10,
    "Dy": 17,
    "Er": 17,
    "Eu": 17,
    "F": 4,
    "Fe": 10,
    "Ga": 9,
    "Gd": 17,
    "Ge": 9,
    "H": 1,
    "He": 1,
    "Hf": 10,
    "Hg": 13,
    "Ho": 17,
    "I": 4,
    "In": 9,
    "Ir": 10,
    "K": 8,
    "Kr": 4,
    "La": 17,
    "Li": 5,
    "Lu": 17,
    "Mg": 4,
    "Mn": 10,
    "Mo": 10,
    "N": 4,
    "Na": 8,
    "Nb": 10,
    "Nd": 17,
    "Ne": 4,
    "Ni": 10,
    "O": 4,
    "Os": 10,
    "P": 4,
    "Pb": 9,
    "Pd": 10,
    "Pm": 17,
    "Po": 9,
    "Pr": 17,
    "Pt": 10,
    "Rb": 8,
    "Re": 10,
    "Rh": 10,
    "Rn": 9,
    "Ru": 10,
    "S": 4,
    "Sb": 9,
    "Sc": 10,
    "Se": 4,
    "Si": 4,
    "Sm": 17,
    "Sn": 9,
    "Sr": 8,
    "Ta": 10,
    "Tb": 17,
    "Tc": 10,
    "Te": 4,
    "Ti": 10,
    "Tl": 9,
    "Tm": 17,
    "V": 10,
    "W": 10,
    "Xe": 9,
    "Y": 10,
    "Yb": 17,
    "Zn": 13,
    "Zr": 10
}
