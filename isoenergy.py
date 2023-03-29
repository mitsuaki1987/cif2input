#!/usr/bin/python3
import os
import pymatgen
from pymatgen.core.periodic_table import get_el_sp
from sssp import psen_dict
import subprocess
import numpy
import sys


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
        structure = pymatgen.core.Structure.from_file(input_file)
        atom = [str(get_el_sp(iat)) for iat in structure.atomic_numbers]
        #
        # Isolated atoms energy
        #
        iso_energy = 0.0
        for iat in atom:
            iso_energy += psen_dict[iat]
        #
        with open(prefix + "_iso.dat", 'w') as f:
            print(iso_energy, file=f)


main()
