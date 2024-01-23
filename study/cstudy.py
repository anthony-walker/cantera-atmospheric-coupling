import os
import re
import sys
import csv
import yaml
import h5py
import numpy
import random
import pandas as pd
import matplotlib.pyplot as plt
from cac.combustor import combustor_atm_sim
from multiprocessing import Pool

cbcolors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']

# change font of plots
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
plt.rcParams["font.size"] = 14

def wrapped(vals):
    outdir = os.path.join(os.path.dirname(__file__), "combustor")
    equiv_ratio, farnesane = vals
    print(f"Running {equiv_ratio}, {farnesane}...")
    log_name = f"combustor-{equiv_ratio:.1f}-{farnesane:.2f}.log"
    log_name = os.path.join("combustor", log_name)
    log_file = open(log_name, "w")
    ostd = sys.stdout
    sys.stdout = log_file
    lf_closed = False
    try:
        combustor_atm_sim(equiv_ratio=equiv_ratio, farnesane=farnesane, outdir=outdir)
    except Exception as e:
        lf_closed = True
        log_file.close()
        sys.stdout = ostd
        print(f"FAILED CASE: {equiv_ratio}, {farnesane}...")
        case_files = filter(lambda x: f"{equiv_ratio:.1f}-{farnesane:.2f}" in x and "log" not in x, os.listdir("combustor"))
        for f in case_files:
            os.remove(os.path.join("combustor", f))
    # close if needed
    if not lf_closed:
        log_file.close()
        sys.stdout = ostd


def run_combustion_study(run_missing=True):
    equiv_ratios = numpy.arange(0.5, 3.0, 0.25)
    farnesane_percents = numpy.arange(0, 0.21, 0.01)
    # get all combinations
    vals_list = []
    for eq in equiv_ratios:
        for xf in farnesane_percents:
            vals_list.append((eq, xf))
    # filter to only missing
    if run_missing:
        all_files = filter(lambda x: "thermo-states" in x, os.listdir("combustor"))
        exists = [re.findall(r"\d+[.]\d+", f) for f in all_files]
        exists.sort()
        keep = []
        for eq, fs in vals_list:
            case = [f"{eq:.1f}", f"{fs:.2f}"]
            if case not in exists:
                keep.append((eq, fs))
        vals_list = keep
    # random.shuffle(vals_list) # random shuffle of order
    # run multi-processing pool
    # equiv_ratio, nsteps, farnesane, outdir
    with Pool(numpy.amin([os.cpu_count(), len(vals_list)])) as p:
        p.map(wrapped, vals_list)

def make_species_contour(unformatted, sp, index=-1):
    equiv_ratios = numpy.arange(0.5, 3.0, 0.25)
    farnesane_percents = numpy.arange(0, 0.21, 0.01)
    masses = numpy.zeros((len(farnesane_percents), len(equiv_ratios)))
    for i, eq in enumerate(equiv_ratios):
        for j, fs in enumerate(farnesane_percents):
            if os.path.exists(os.path.join("combustor", unformatted.format(eq, fs))):
                short_data = pd.read_csv(os.path.join("combustor", unformatted.format(eq, fs)), delimiter=",", engine="python", usecols=["mass", "time", f'Y_{sp}'])
                masses[j][i] = short_data[f'Y_{sp}'].iloc[index] * short_data["mass"].iloc[index]
            else:
                masses[j][i] = numpy.NAN
    # make contour
    X, Y = numpy.meshgrid(equiv_ratios, farnesane_percents)
    plt.contourf(X, Y, masses, cmap='inferno')
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Farnesane Percentage')
    plt.colorbar(label='Mass [kg]')
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.savefig(f"figures/{sp}-contour.pdf")

def short_term_states_plots(species):
    short_data = pd.read_csv(os.path.join("combustor", 'short-term-states-1.0-0.00.csv'), delimiter=",", engine="python")
    mass = short_data["mass"]
    time = short_data["time"]
    Y_sp = short_data[f"Y_{species}"]
    fig, ax = plt.subplots()
    # ax.plot(time, mass * Y_tol)
    ax.plot(time, mass * Y_sp)
    # ax.plot(time, mass)
    plt.show()

if __name__ == "__main__":
    # run_combustion_study()
    make_species_contour("short-term-states-{:.1f}-{:.2f}.csv", "TOLUENE")
    make_species_contour("long-term-states-{:.1f}-{:.2f}.csv", "TOLUENE")

