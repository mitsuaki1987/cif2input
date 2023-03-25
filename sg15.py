import numpy


pseudo_dict = {
    "Ag": "Ag_ONCV_PBE-1.0.upf",
    "Al": "Al_ONCV_PBE-1.0.upf",
    "Ar": "Ar_ONCV_PBE-1.1.upf",
    "As": "As_ONCV_PBE-1.1.upf",
    "Au": "Au_ONCV_PBE-1.0.upf",
    "B": "B_ONCV_PBE-1.0.upf",
    "Ba": "Ba_ONCV_PBE-1.0.upf",
    "Be": "Be_ONCV_PBE-1.0.upf",
    "Bi": "Bi_ONCV_PBE-1.0.upf",
    "Br": "Br_ONCV_PBE-1.0.upf",
    "C": "C_ONCV_PBE-1.0.upf",
    "Ca": "Ca_ONCV_PBE-1.0.upf",
    "Cd": "Cd_ONCV_PBE-1.0.upf",
    "Cl": "Cl_ONCV_PBE-1.1.upf",
    "Co": "Co_ONCV_PBE-1.0.upf",
    "Cr": "Cr_ONCV_PBE-1.0.upf",
    "Cs": "Cs_ONCV_PBE-1.0.upf",
    "Cu": "Cu_ONCV_PBE-1.0.upf",
    "F": "F_ONCV_PBE-1.0.upf",
    "Fe": "Fe_ONCV_PBE-1.0.upf",
    "Ga": "Ga_ONCV_PBE-1.0.upf",
    "Ge": "Ge_ONCV_PBE-1.0.upf",
    "H": "H_ONCV_PBE-1.0.upf",
    "He": "He_ONCV_PBE-1.0.upf",
    "Hf": "Hf_ONCV_PBE-1.0.upf",
    "Hg": "Hg_ONCV_PBE-1.0.upf",
    "I": "I_ONCV_PBE-1.1.upf",
    "In": "In_ONCV_PBE-1.1.upf",
    "Ir": "Ir_ONCV_PBE-1.0.upf",
    "K": "K_ONCV_PBE-1.0.upf",
    "Kr": "Kr_ONCV_PBE-1.0.upf",
    "La": "La_ONCV_PBE-1.0.upf",
    "Li": "Li_ONCV_PBE-1.0.upf",
    "Mg": "Mg_ONCV_PBE-1.0.upf",
    "Mn": "Mn_ONCV_PBE-1.0.upf",
    "Mo": "Mo_ONCV_PBE-1.0.upf",
    "N": "N_ONCV_PBE-1.0.upf",
    "Na": "Na_ONCV_PBE-1.0.upf",
    "Nb": "Nb_ONCV_PBE-1.0.upf",
    "Ne": "Ne_ONCV_PBE-1.0.upf",
    "Ni": "Ni_ONCV_PBE-1.0.upf",
    "O": "O_ONCV_PBE-1.0.upf",
    "Os": "Os_ONCV_PBE-1.0.upf",
    "P": "P_ONCV_PBE-1.1.upf",
    "Pb": "Pb_ONCV_PBE-1.0.upf",
    "Pd": "Pd_ONCV_PBE-1.0.upf",
    "Pt": "Pt_ONCV_PBE-1.0.upf",
    "Rb": "Rb_ONCV_PBE-1.0.upf",
    "Re": "Re_ONCV_PBE-1.0.upf",
    "Rh": "Rh_ONCV_PBE-1.0.upf",
    "Ru": "Ru_ONCV_PBE-1.0.upf",
    "S": "S_ONCV_PBE-1.0.upf",
    "Sb": "Sb_ONCV_PBE-1.1.upf",
    "Sc": "Sc_ONCV_PBE-1.0.upf",
    "Se": "Se_ONCV_PBE-1.1.upf",
    "Si": "Si_ONCV_PBE-1.0.upf",
    "Sn": "Sn_ONCV_PBE-1.1.upf",
    "Sr": "Sr_ONCV_PBE-1.0.upf",
    "Ta": "Ta_ONCV_PBE-1.0.upf",
    "Tc": "Tc_ONCV_PBE-1.0.upf",
    "Te": "Te_ONCV_PBE-1.0.upf",
    "Ti": "Ti_ONCV_PBE-1.0.upf",
    "Tl": "Tl_ONCV_PBE-1.0.upf",
    "V": "V_ONCV_PBE-1.0.upf",
    "W": "W_ONCV_PBE-1.0.upf",
    "Xe": "Xe_ONCV_PBE-1.1.upf",
    "Y": "Y_ONCV_PBE-1.0.upf",
    "Zn": "Zn_ONCV_PBE-1.0.upf",
    "Zr": "Zr_ONCV_PBE-1.0.upf"
}

ecutwfc_dict = {
    "Ag": 50.0,
    "Al": 65.0,
    "Ar": 60.0,
    "As": 35.0,
    "Au": 45.0,
    "B": 40.0,
    "Ba": 30.0,
    "Be": 65.0,
    "Bi": 120.0,
    "Br": 80.0,
    "C": 65.0,
    "Ca": 120.0,
    "Cd": 45.0,
    "Cl": 75.0,
    "Co": 60.0,
    "Cr": 50.0,
    "Cs": 75.0,
    "Cu": 90.0,
    "F": 75.0,
    "Fe": 120.0,
    "Ga": 150.0,
    "Ge": 70.0,
    "H": 70.0,
    "He": 50.0,
    "Hf": 55.0,
    "Hg": 50.0,
    "I": 80.0,
    "In": 65.0,
    "Ir": 40.0,
    "K": 120.0,
    "Kr": 45.0,
    "La": 60.0,
    "Li": 60.0,
    "Mg": 65.0,
    "Mn": 60.0,
    "Mo": 35.0,
    "N": 75.0,
    "Na": 90.0,
    "Nb": 90.0,
    "Ne": 200.0,
    "Ni": 65.0,
    "O": 120.0,
    "Os": 55.0,
    "P": 55.0,
    "Pb": 35.0,
    "Pd": 45.0,
    "Pt": 60.0,
    "Rb": 30.0,
    "Re": 60.0,
    "Rh": 35.0,
    "Ru": 35.0,
    "S": 55.0,
    "Sb": 60.0,
    "Sc": 45.0,
    "Se": 30.0,
    "Si": 45.0,
    "Sn": 50.0,
    "Sr": 30.0,
    "Ta": 50.0,
    "Tc": 30.0,
    "Te": 70.0,
    "Ti": 50.0,
    "Tl": 55.0,
    "V": 100.0,
    "W": 50.0,
    "Xe": 60.0,
    "Y": 40.0,
    "Zn": 90.0,
    "Zr": 50.0
}

ecutrho_dict = {
    "Ag": 200.0,
    "Al": 260.0,
    "Ar": 240.0,
    "As": 140.0,
    "Au": 180.0,
    "B": 160.0,
    "Ba": 120.0,
    "Be": 260.0,
    "Bi": 480.0,
    "Br": 320.0,
    "C": 260.0,
    "Ca": 480.0,
    "Cd": 180.0,
    "Cl": 300.0,
    "Co": 240.0,
    "Cr": 200.0,
    "Cs": 300.0,
    "Cu": 360.0,
    "F": 300.0,
    "Fe": 480.0,
    "Ga": 600.0,
    "Ge": 280.0,
    "H": 280.0,
    "He": 200.0,
    "Hf": 220.0,
    "Hg": 200.0,
    "I": 320.0,
    "In": 260.0,
    "Ir": 160.0,
    "K": 480.0,
    "Kr": 180.0,
    "La": 240.0,
    "Li": 240.0,
    "Mg": 260.0,
    "Mn": 240.0,
    "Mo": 140.0,
    "N": 300.0,
    "Na": 360.0,
    "Nb": 360.0,
    "Ne": 800.0,
    "Ni": 260.0,
    "O": 480.0,
    "Os": 220.0,
    "P": 220.0,
    "Pb": 140.0,
    "Pd": 180.0,
    "Pt": 240.0,
    "Rb": 120.0,
    "Re": 240.0,
    "Rh": 140.0,
    "Ru": 140.0,
    "S": 220.0,
    "Sb": 240.0,
    "Sc": 180.0,
    "Se": 120.0,
    "Si": 180.0,
    "Sn": 200.0,
    "Sr": 120.0,
    "Ta": 200.0,
    "Tc": 120.0,
    "Te": 280.0,
    "Ti": 200.0,
    "Tl": 220.0,
    "V": 400.0,
    "W": 200.0,
    "Xe": 240.0,
    "Y": 160.0,
    "Zn": 360.0,
    "Zr": 200.0
}

valence_dict = {
    "Ag": 19,
    "Al": 11,
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
    "Cd": 20,
    "Cl": 7,
    "Co": 17,
    "Cr": 14,
    "Cs": 9,
    "Cu": 19,
    "F": 7,
    "Fe": 16,
    "Ga": 13,
    "Ge": 14,
    "H": 1,
    "He": 2,
    "Hf": 26,
    "Hg": 20,
    "I": 17,
    "In": 13,
    "Ir": 17,
    "K": 9,
    "Kr": 8,
    "La": 11,
    "Li": 3,
    "Mg": 10,
    "Mn": 15,
    "Mo": 14,
    "N": 5,
    "Na": 9,
    "Nb": 13,
    "Ne": 8,
    "Ni": 18,
    "O": 6,
    "Os": 16,
    "P": 5,
    "Pb": 14,
    "Pd": 18,
    "Pt": 18,
    "Rb": 9,
    "Re": 15,
    "Rh": 17,
    "Ru": 16,
    "S": 6,
    "Sb": 15,
    "Sc": 11,
    "Se": 6,
    "Si": 4,
    "Sn": 14,
    "Sr": 10,
    "Ta": 27,
    "Tc": 15,
    "Te": 16,
    "Ti": 12,
    "Tl": 13,
    "V": 13,
    "W": 28,
    "Xe": 18,
    "Y": 11,
    "Zn": 20,
    "Zr": 12
}

atomwfc_dict = {
    "Ag": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Al": numpy.array([["2S", 1], ["2P", 3], ["3S", 1], ["3P", 3]]),
    "Ar": numpy.array([["3S", 1], ["3P", 3]]),
    "As": numpy.array([["4S", 1], ["4P", 3]]),
    "Au": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "B": numpy.array([["2S", 1], ["2P", 3]]),
    "Ba": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1]]),
    "Be": numpy.array([["1S", 1], ["2S", 1]]),
    "Bi": numpy.array([["5D", 5], ["6S", 1], ["6P", 3]]),
    "Br": numpy.array([["4S", 1], ["4P", 3]]),
    "C": numpy.array([["2S", 1], ["2P", 3]]),
    "Ca": numpy.array([["3S", 1], ["3P", 3], ["4S", 1]]),
    "Cd": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Cl": numpy.array([["3S", 1], ["3P", 3]]),
    "Co": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Cr": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Cs": numpy.array([["5S", 1], ["5P", 3], ["6S", 1]]),
    "Cu": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "F": numpy.array([["2S", 1], ["2P", 3]]),
    "Fe": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Ga": numpy.array([["3D", 5], ["4S", 1], ["4P", 3]]),
    "Ge": numpy.array([["3D", 5], ["4S", 1], ["4P", 3]]),
    "H": numpy.array([["1S", 1]]),
    "He": numpy.array([["1S", 1]]),
    "Hf": numpy.array([["4F", 7], ["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "Hg": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "I": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "In": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Ir": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "K": numpy.array([["3S", 1], ["3P", 3], ["4S", 1]]),
    "Kr": numpy.array([["4S", 1], ["4P", 3]]),
    "La": numpy.array([["5S", 1], ["5P", 3], ["5D", 5], ["6S", 1]]),
    "Li": numpy.array([["1S", 1], ["2S", 1]]),
    "Mg": numpy.array([["2S", 1], ["2P", 3], ["3S", 1]]),
    "Mn": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Mo": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "N": numpy.array([["2S", 1], ["2P", 3]]),
    "Na": numpy.array([["2S", 1], ["2P", 3], ["3S", 1]]),
    "Nb": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Ne": numpy.array([["2S", 1], ["2P", 3]]),
    "Ni": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "O": numpy.array([["2S", 1], ["2P", 3]]),
    "Os": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "P": numpy.array([["3S", 1], ["3P", 3]]),
    "Pb": numpy.array([["5D", 5], ["6S", 1], ["6P", 3]]),
    "Pd": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Pt": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "Rb": numpy.array([["4S", 1], ["4P", 3], ["5S", 1]]),
    "Re": numpy.array([["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "Rh": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Ru": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "S": numpy.array([["3S", 1], ["3P", 3]]),
    "Sb": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Sc": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Se": numpy.array([["4S", 1], ["4P", 3]]),
    "Si": numpy.array([["3S", 1], ["3P", 3]]),
    "Sn": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Sr": numpy.array([["4S", 1], ["4P", 3], ["5S", 1]]),
    "Ta": numpy.array([["4F", 7], ["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "Tc": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Te": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Ti": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Tl": numpy.array([["5D", 5], ["6S", 1], ["6P", 3]]),
    "V": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "W": numpy.array([["4F", 7], ["5S", 1], ["5P", 3], ["6S", 1], ["5D", 5]]),
    "Xe": numpy.array([["4D", 5], ["5S", 1], ["5P", 3]]),
    "Y": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]]),
    "Zn": numpy.array([["3S", 1], ["3P", 3], ["4S", 1], ["3D", 5]]),
    "Zr": numpy.array([["4S", 1], ["4P", 3], ["5S", 1], ["4D", 5]])
}

psen_dict = {
    "Ag": -2.87458497570e+02,
    "Al": -1.38158630160e+02,
    "Ar": -4.23888901800e+01,
    "As": -1.23803319800e+01,
    "Au": -2.74500089190e+02,
    "B": -5.19895176000e+00,
    "Ba": -5.13112242500e+01,
    "Be": -2.63332978400e+01,
    "Bi": -1.35667183940e+02,
    "Br": -2.68737806900e+01,
    "C": -1.07608311900e+01,
    "Ca": -7.35908384900e+01,
    "Cd": -3.10809138440e+02,
    "Cl": -3.00778375800e+01,
    "Co": -2.70558259780e+02,
    "Cr": -1.66964014250e+02,
    "Cs": -4.03893055200e+01,
    "Cu": -3.64938847510e+02,
    "F": -4.74923082700e+01,
    "Fe": -2.36173028060e+02,
    "Ga": -1.31110149120e+02,
    "Ge": -1.41282620400e+02,
    "H": -9.31547810000e-01,
    "He": -5.58120348000e+00,
    "Hf": -4.46943719610e+02,
    "Hg": -3.10268913870e+02,
    "I": -1.98592398810e+02,
    "In": -1.04422266250e+02,
    "Ir": -2.11891726710e+02,
    "K": -5.67680287900e+01,
    "Kr": -3.71776222000e+01,
    "La": -6.35722511700e+01,
    "Li": -1.41218982800e+01,
    "Mg": -1.08540968200e+02,
    "Mn": -1.98506841210e+02,
    "Mo": -1.36441981150e+02,
    "N": -1.93311748300e+01,
    "Na": -8.51189606200e+01,
    "Nb": -1.13494989010e+02,
    "Ne": -6.13600994400e+01,
    "Ni": -3.07951024790e+02,
    "O": -3.14735372000e+01,
    "Os": -1.83322536990e+02,
    "P": -1.28414963200e+01,
    "Pb": -1.12671379520e+02,
    "Pd": -2.53104481200e+02,
    "Pt": -2.42207146750e+02,
    "Rb": -4.88737685900e+01,
    "Re": -1.59067573290e+02,
    "Rh": -2.16669801870e+02,
    "Ru": -1.88069238110e+02,
    "S": -2.02254592100e+01,
    "Sb": -1.40809074850e+02,
    "Sc": -9.15691204500e+01,
    "Se": -1.86879641300e+01,
    "Si": -7.52258959000e+00,
    "Sn": -1.20405401260e+02,
    "Sr": -6.21401700800e+01,
    "Ta": -4.96006509420e+02,
    "Tc": -1.60575133670e+02,
    "Te": -1.75629608210e+02,
    "Ti": -1.14938194020e+02,
    "Tl": -9.57793795600e+01,
    "V": -1.41694200160e+02,
    "W": -5.43194215200e+02,
    "Xe": -2.31787053180e+02,
    "Y": -7.68704464600e+01,
    "Zn": -4.00242728200e+02,
    "Zr": -9.43829629500e+01
}

band_dict = {
    "Ag": 10,
    "Al": 8,
    "Ar": 4,
    "As": 4,
    "Au": 10,
    "B": 4,
    "Ba": 10,
    "Be": 2,
    "Bi": 9,
    "Br": 4,
    "C": 4,
    "Ca": 5,
    "Cd": 10,
    "Cl": 4,
    "Co": 10,
    "Cr": 10,
    "Cs": 5,
    "Cu": 10,
    "F": 4,
    "Fe": 10,
    "Ga": 9,
    "Ge": 9,
    "H": 1,
    "He": 1,
    "Hf": 17,
    "Hg": 10,
    "I": 9,
    "In": 9,
    "Ir": 10,
    "K": 5,
    "Kr": 4,
    "La": 10,
    "Li": 2,
    "Mg": 5,
    "Mn": 10,
    "Mo": 10,
    "N": 4,
    "Na": 5,
    "Nb": 10,
    "Ne": 4,
    "Ni": 10,
    "O": 4,
    "Os": 10,
    "P": 4,
    "Pb": 9,
    "Pd": 10,
    "Pt": 10,
    "Rb": 4,
    "Re": 10,
    "Rh": 10,
    "Ru": 10,
    "S": 4,
    "Sb": 9,
    "Sc": 10,
    "Se": 4,
    "Si": 4,
    "Sn": 9,
    "Sr": 5,
    "Ta": 17,
    "Tc": 10,
    "Te": 9,
    "Ti": 10,
    "Tl": 9,
    "V": 10,
    "W": 17,
    "Xe": 9,
    "Y": 10,
    "Zn": 10,
    "Zr": 10
}

core_dict = {
    "Ag": 4,
    "Al": 4,
    "Ar": 0,
    "As": 0,
    "Au": 4,
    "B": 0,
    "Ba": 4,
    "Be": 1,
    "Bi": 5,
    "Br": 0,
    "C": 0,
    "Ca": 4,
    "Cd": 4,
    "Cl": 0,
    "Co": 4,
    "Cr": 4,
    "Cs": 4,
    "Cu": 4,
    "F": 0,
    "Fe": 4,
    "Ga": 5,
    "Ge": 5,
    "H": 0,
    "He": 0,
    "Hf": 11,
    "Hg": 4,
    "I": 5,
    "In": 5,
    "Ir": 4,
    "K": 4,
    "Kr": 0,
    "La": 4,
    "Li": 1,
    "Mg": 4,
    "Mn": 4,
    "Mo": 4,
    "N": 0,
    "Na": 4,
    "Nb": 4,
    "Ne": 0,
    "Ni": 4,
    "O": 0,
    "Os": 4,
    "P": 0,
    "Pb": 5,
    "Pd": 4,
    "Pt": 4,
    "Rb": 4,
    "Re": 4,
    "Rh": 4,
    "Ru": 4,
    "S": 0,
    "Sb": 5,
    "Sc": 4,
    "Se": 0,
    "Si": 0,
    "Sn": 5,
    "Sr": 4,
    "Ta": 11,
    "Tc": 4,
    "Te": 5,
    "Ti": 4,
    "Tl": 5,
    "V": 4,
    "W": 11,
    "Xe": 5,
    "Y": 4,
    "Zn": 4,
    "Zr": 4
}

wan_dict = {
    "Ag": ["s", "d"],
    "Al": ["s", "p"],
    "Ar": ["s", "p"],
    "As": ["s", "p"],
    "Au": ["s", "d"],
    "B": ["s", "p"],
    "Ba": ["s", "d"],
    "Be": ["s"],
    "Bi": ["s", "p"],
    "Br": ["s", "p"],
    "C": ["s", "p"],
    "Ca": ["s"],
    "Cd": ["s", "d"],
    "Cl": ["s", "p"],
    "Co": ["s", "d"],
    "Cr": ["s", "d"],
    "Cs": ["s"],
    "Cu": ["s", "d"],
    "F": ["s", "p"],
    "Fe": ["s", "d"],
    "Ga": ["s", "p"],
    "Ge": ["s", "p"],
    "H": ["s"],
    "He": ["s"],
    "Hf": ["s", "d"],
    "Hg": ["s", "d"],
    "I": ["s", "p"],
    "In": ["s", "p"],
    "Ir": ["s", "d"],
    "K": ["s"],
    "Kr": ["s", "p"],
    "La": ["s", "d"],
    "Li": ["s"],
    "Mg": ["s"],
    "Mn": ["s", "d"],
    "Mo": ["s", "d"],
    "N": ["s", "p"],
    "Na": ["s"],
    "Nb": ["s", "d"],
    "Ne": ["s", "p"],
    "Ni": ["s", "d"],
    "O": ["s", "p"],
    "Os": ["s", "d"],
    "P": ["s", "p"],
    "Pb": ["s", "p"],
    "Pd": ["s", "d"],
    "Pt": ["s", "d"],
    "Rb": ["s"],
    "Re": ["s", "d"],
    "Rh": ["s", "d"],
    "Ru": ["s", "d"],
    "S": ["s", "p"],
    "Sb": ["s", "p"],
    "Sc": ["s", "d"],
    "Se": ["s", "p"],
    "Si": ["s", "p"],
    "Sn": ["s", "p"],
    "Sr": ["s"],
    "Ta": ["s", "d"],
    "Tc": ["s", "d"],
    "Te": ["s", "p"],
    "Ti": ["s", "d"],
    "Tl": ["s", "p"],
    "V": ["s", "d"],
    "W": ["s", "d"],
    "Xe": ["s", "p"],
    "Y": ["s", "d"],
    "Zn": ["s", "d"],
    "Zr": ["s", "d"]
}
