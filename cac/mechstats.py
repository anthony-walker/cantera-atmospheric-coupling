import os
import sys
import cantera as ct
import warnings
from cac.constants import DATA_DIR


def all_mechanisms():
    simple = os.path.join(DATA_DIR, "sulfur", "sulfur.yaml")
    print_mech_stats(simple, "gas")

    simple = os.path.join(DATA_DIR, "farnesane", "farnesane.yaml")
    print_mech_stats(simple, "gas")

    sc = os.path.join(DATA_DIR, "combustor-min.yaml")
    print_mech_stats(sc, "combustor")

    sc = os.path.join(DATA_DIR, "combustor-min.yaml")
    print_mech_stats(sc, "atmosphere")


def print_mech_stats(fname, phasename):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gas = ct.Solution(fname, phasename, transport_model=None)
        print(os.path.basename(fname))
        print(phasename)
        print(gas.n_species)
        print(gas.n_reactions)
        print()
