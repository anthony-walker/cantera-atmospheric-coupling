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
from cac.constants import COLORS

MASK_VALUE = 1e-40
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

def make_species_contour(sp, index=-1, path="fullmcm", name="", mode="hdf5", scale="short"):
    equiv_ratios = numpy.arange(0.7, 2.6, 0.1)
    farnesane_percents = numpy.arange(0, 0.21, 0.01)
    unformatted = f"{scale}-term-states-"+ "{:.1f}-{:.2f}" + f".{mode}"
    masses = numpy.zeros((len(farnesane_percents), len(equiv_ratios)))
    for i, eq in enumerate(equiv_ratios):
        for j, fs in enumerate(farnesane_percents):
            if os.path.exists(os.path.join(path, unformatted.format(eq, fs))):
                if mode == "csv":
                    short_data = pd.read_csv(os.path.join(path, unformatted.format(eq, fs)), delimiter=",", engine="python", usecols=["mass", "time", f'Y_{sp}'])
                    masses[j][i] = short_data[f'Y_{sp}'].iloc[index] * short_data['mass'].iloc[index]
                elif mode == "hdf5":
                    short_data = h5py.File(os.path.join(path, unformatted.format(eq, fs)), "r")
                    masses[j][i] = short_data[f'Y_{sp}'][index] * short_data['mass'][index]
            else:
                masses[j][i] = numpy.NAN
    # apply mask to masses
    mask = numpy.abs(masses) <= MASK_VALUE
    for i in range(len(equiv_ratios)):
        for j in range(len(farnesane_percents)):
            masses[j][i] = masses[j][i] if not mask[j][i] else 0
    # make contour
    X, Y = numpy.meshgrid(equiv_ratios, farnesane_percents)
    plt.contourf(X, Y, masses, cmap='viridis')
    plt.xticks(numpy.linspace(0.7, 2.5, 5))
    plt.yticks(numpy.linspace(0, 0.2, 5))
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Farnesane Fraction')
    plt.colorbar(label='Mass [kg]')
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/{sp.lower()}-{scale}-{index}-contour.pdf", bbox_inches='tight')
    plt.close()


def states_plots(species=[], path="fullmcm", mode="csv", scale="long", eq=1.0, fp=0.10, name=None, normalized=False, limit=None):
    if mode == "csv":
        short_data = pd.read_csv(os.path.join(path, f'{scale}-term-states-{eq:.1f}-{fp:.2f}.csv'), delimiter=",", engine="python")
        mass = short_data["mass"]
        time = short_data["time"]
    elif mode == "hdf5":
        short_data = h5py.File(os.path.join(path, f'{scale}-term-states-{eq:.1f}-{fp:.2f}.hdf5'), "r")
        mass = numpy.array(short_data["mass"])
        time = numpy.array(short_data["time"])
    # adjust time scale
    if scale == "long":
        time = time // 24 // 3600 # days
    if species is None:
        for k in short_data.keys():
            data = numpy.array(short_data[k])
            mask = data <= MASK_VALUE
            for i in range(len(data)):
                data[i] = data[i] if not mask[i] else 0

            if limit is None or numpy.amax(data) > limit:
                if normalized:
                    data = data / numpy.amax(data)
                fig, ax = plt.subplots()
                ax.plot(time, data)
                ax.set_title(k)
                ax.set_yticks(numpy.linspace(0, numpy.amax(data) * 1.1, 5))
                ax.set_ylabel("Mass [kg]")
                if scale == "long":
                    ax.set_xlabel("Time [days]")
                else:
                    ax.set_xlabel("Time [s]")
                plt.tight_layout()
                plt.savefig(f"figures/{k}-{scale}.pdf")
                plt.close()
    else:
        if isinstance(species, str):
            species = [species]
        fig, ax = plt.subplots()
        for j, sp in enumerate(species):
            data = numpy.array(short_data[f"Y_{sp}"])
            mask = numpy.abs(data) <= MASK_VALUE
            for i in range(len(data)):
                data[i] = data[i] if not mask[i] else 0
            if normalized:
                data = data / numpy.amax(data)
            lb = sp.lower() if sp != "C5H8" else "isoprene"
            ax.plot(time, data, color=COLORS[j], label=lb, linewidth=3)
        fig.legend(loc='upper center', ncols=len(species)//2, bbox_to_anchor=(0.5, 1.125))
        if scale == "short":
            plt.xlabel('Time [s]')
        else:
            plt.xlabel("Time [days]")
        if normalized:
            plt.ylabel('Normalized Mass')
        else:
            plt.ylabel('Mass [kg]')
        plt.grid(True)
        plt.tight_layout()
        if name is not None:
            plt.savefig(f"figures/{name}-{scale}-{eq:.1f}-{fp:.2f}.pdf", bbox_inches='tight')
        else:
            plt.savefig(f"figures/states-{scale}-{eq:.1f}-{fp:.2f}.pdf", bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    mode = "hdf5"
    path = "maximal"
    all_specs = ["TOLUENE", "PHENOL", "BENZENE", "C5H8"]
    states_plots(None, path=path, mode=mode, scale="long", eq=1.5, fp=0.08, limit=1e-24)
    states_plots(all_specs, path=path, mode=mode, scale="short", eq=0.7, fp=0.08, normalized=True)
    states_plots(all_specs, path=path, mode=mode, scale="short", eq=1.5, fp=0.08, normalized=True)

    for name in all_specs: #
        make_species_contour(name, index=0, path=path, mode=mode)
        # make_species_contour(name, index=0, path=path, mode=mode, scale="long")
    peroxys = ["BZBIPERO2", "TLBIPERO2"]
    states_plots(peroxys, path=path, mode=mode, scale="long", eq=1.5, fp=0.08, name=f"peroxys-states", normalized=True)
    for ro2 in peroxys:
        make_species_contour(ro2, path=path, scale="long", mode=mode)
        states_plots(ro2, path=path, mode=mode, scale="long", eq=1.5, fp=0.08, name=f"{ro2.lower()}-states")
