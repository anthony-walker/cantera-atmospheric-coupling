import os
import sys
import cantera as ct
import warnings
from cac.constants import DATA_DIR


def all_mechanisms():
    simple = os.path.join(DATA_DIR, "simple", "simple.yaml")
    print_mech_stats(simple, "atmosphere")

    simple = os.path.join(DATA_DIR, "atmosphere", "atmosphere.yaml")
    print_mech_stats(simple, "atmosphere")

    simple = os.path.join(DATA_DIR, "sulfur", "sulfur.yaml")
    print_mech_stats(simple, "gas")

    simple = os.path.join(DATA_DIR, "creck", "jetfuel.yaml")
    print_mech_stats(simple, "gas")

    simple = os.path.join(DATA_DIR, "farnesene", "farnesene-313-2147.yaml")
    print_mech_stats(simple, "gas")

    sc = os.path.join(DATA_DIR, "scombustor.yaml")
    print_mech_stats(sc, "gas")


def print_mech_stats(fname, phasename):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gas = ct.Solution(fname, phasename, transport_model=None)
        print(os.path.basename(fname))
        print(phasename)
        print(gas.n_species)
        print(gas.n_reactions)
        print()
