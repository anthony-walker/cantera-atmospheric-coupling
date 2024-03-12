import os
import re
import sys
import csv
import copy
import yaml
import math
import h5py
import numpy
import random
import warnings
import cantera as ct
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from cac.verification import polimi_model
from cac.combustor import combustor_atm_sim, multizone_combustor
from multiprocessing import Pool
from cac.constants import COLORS, DATA_DIR, MARKERS, HATCHES
from cac.reactors import PlumeSolution
from openap import FuelFlow, Emission
from matplotlib.lines import Line2D

MASK_VALUE = 1e-20
cbcolors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']

# change font of plots
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
plt.rcParams["font.size"] = 16

def make_temperature_contour(index=-1, path="fullmcm", mode="hdf5", scale="short"):
    equiv_ratios = numpy.arange(0.2, 1.1, 0.1)
    farnesane_percents = numpy.arange(0, 0.21, 0.01)
    unformatted = f"{scale}-term-states-"+ "{:.1f}-{:.2f}" + f".{mode}"
    yf = "thermo-states-{:.1f}-{:.2f}.yaml"
    temps = numpy.zeros((len(farnesane_percents), len(equiv_ratios)))
    for i, eq in enumerate(equiv_ratios):
        for j, fs in enumerate(farnesane_percents):
            if os.path.exists(os.path.join(path, unformatted.format(eq, fs))) and mode != "yaml":
                if mode == "csv":
                    short_data = pd.read_csv(os.path.join(path, unformatted.format(eq, fs)), delimiter=",", engine="python", usecols=[f'T'])
                    temps[j][i] = short_data[f'T'].iloc[index]
                elif mode == "hdf5":
                    short_data = h5py.File(os.path.join(path, unformatted.format(eq, fs)), "r")
                    temps[j][i] = short_data[f'T'][index]
            elif os.path.exists(os.path.join(path, yf.format(eq, fs))):
                with open(os.path.join(path, yf.format(eq, fs)), "r") as f:
                    yaml_data = yaml.load(f, Loader=yaml.SafeLoader)
                temps[j][i] = float(yaml_data["T"][index])
            else:
                temps[j][i] = numpy.NAN
    # apply mask to temps
    mask = numpy.abs(temps) <= MASK_VALUE
    for i in range(len(equiv_ratios)):
        for j in range(len(farnesane_percents)):
            temps[j][i] = temps[j][i] if not mask[j][i] else 0
    # make contour
    X, Y = numpy.meshgrid(equiv_ratios, farnesane_percents)
    plt.contourf(X, Y, temps, cmap='viridis')
    plt.xticks(numpy.linspace(0.2, 1.0, 5))
    plt.yticks(numpy.linspace(0, 0.2, 5))
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Farnesane Fraction')
    plt.colorbar(label='Temperature [K]')
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/{path}-temperature-{scale}-{index}-contour.pdf", bbox_inches='tight')
    plt.close()

def make_species_contour(sp, index=-1, path="fullmcm", name="", mode="hdf5", scale="long-term", ekey="mdot_exhaust"):
    equiv_ratios = numpy.arange(0.2, 1.1, 0.1)
    farnesane_percents = numpy.arange(0, 0.21, 0.01)
    unformatted = f"{scale}-states-"+ "{:.1f}-{:.2f}" + f".{mode}"
    masses = numpy.zeros((len(farnesane_percents), len(equiv_ratios)))
    for i, eq in enumerate(equiv_ratios):
        for j, fs in enumerate(farnesane_percents):
            if os.path.exists(os.path.join(path, unformatted.format(eq, fs))):
                # load yaml file for rate
                with open(os.path.join(path, f"thermo-states-{eq:.1f}-{fs:.2f}.yaml"), "r") as tsf:
                    yaml_data = yaml.load(tsf, Loader=yaml.SafeLoader)

                if mode == "csv":
                    short_data = pd.read_csv(os.path.join(path, unformatted.format(eq, fs)), delimiter=",", engine="python", usecols=["mass", "time", f'Y_{sp}'])
                    masses[j][i] = short_data[f'Y_{sp}'].iloc[index] * yaml_data[ekey] / short_data[f'mass'].iloc[index] * 1e6
                elif mode == "hdf5":
                    short_data = h5py.File(os.path.join(path, unformatted.format(eq, fs)), "r")
                    masses[j][i] = short_data[f'Y_{sp}'][index] * yaml_data[ekey] / short_data[f'mass'][index] * 1e6
            else:
                masses[j][i] = numpy.NAN
    # apply mask to masses
    mask = numpy.abs(masses) <= MASK_VALUE
    mask = numpy.logical_or(numpy.isnan(masses), mask)
    for i in range(len(equiv_ratios)):
        for j in range(len(farnesane_percents)):
            masses[j][i] = masses[j][i] if not mask[j][i] else MASK_VALUE/10
    # make contour
    X, Y = numpy.meshgrid(equiv_ratios, farnesane_percents)
    fig, ax = plt.subplots()
    fig.set_size_inches(5.5, 8)
    contour = plt.contourf(X, Y, masses, cmap='viridis', norm=LogNorm(vmin=MASK_VALUE/10, vmax=numpy.amax(masses)))
    plt.xticks(numpy.linspace(0.2, 1.0, 5))
    plt.yticks(numpy.linspace(0, 0.2, 5))
    plt.xlabel('Equivalence Ratio', fontsize=20)
    plt.ylabel('Farnesane Fraction', fontsize=20)
    cbar = plt.colorbar(label=r'Production [$\text{ppm} / \text{s}$]')
    cbar.ax.tick_params(labelsize=24)
    cbar.set_label(r'Production [$\text{ppm} / \text{s}$]', fontsize=20)
    # Get the current tick labels
    ticks = [t for t in cbar.get_ticks()]
    ticks = numpy.logspace(numpy.log10(ticks[0]), numpy.log10(ticks[-1]), 6)
    tick_labels = ["$10^{" + f"{int(numpy.log10(i)):d}" + "}$" for i in ticks]
    # Modify only the first tick label
    tick_labels[0] = '$0$'
    # Set the modified tick labels
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(tick_labels)
    # mkdir if it does not exist
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/{path}-{sp.lower()}-{scale}-{index}-contour.pdf", bbox_inches='tight')
    plt.close()


def states_plots(species=[], path="fullmcm", mode="csv", scale="long", eq=1.0, fp=0.10, name=None, normalized=False, limit=MASK_VALUE, pltfcn=plt.plot, ncols=2, markers=[], mf=False, legend_on=True, return_fig=False, fig=None, ax=None, linestyle="-", entr_on=True, CID=0):
    # load data
    with open(os.path.join(path, f"thermo-states-{eq:.1f}-{fp:.2f}.yaml"), "r") as tsf:
        mdot_exhaust = yaml.load(tsf, Loader=yaml.SafeLoader)["mdot_exhaust"]
    if mode == "csv":
        short_data = pd.read_csv(os.path.join(path, f'{scale}-term-states-{eq:.1f}-{fp:.2f}.csv'), delimiter=",", engine="python")
        mass = short_data["mass"]
        time = short_data["time"]
    elif mode == "hdf5":
        short_data = h5py.File(os.path.join(path, f'{scale}-term-states-{eq:.1f}-{fp:.2f}.hdf5'), "r")
        mass = numpy.array(short_data["mass"])
        time = numpy.array(short_data["time"])
    # adjust time scale
    # if scale == "long":
    #     time = time // 24 // 3600 # days
    if species is None:
        for k in short_data.keys():
            data = numpy.array(short_data[k]) * mdot_exhaust / numpy.array(mass) * 1e6
            mask = data <= MASK_VALUE
            for i in range(len(data)):
                data[i] = data[i] if not mask[i] else 0

            if limit is None or numpy.amax(data) > limit:
                if normalized:
                    data = data / numpy.amax(data)
                fig, ax = plt.subplots()
                pltfcn(time, data)
                ax.set_title(k)
                # ax.set_yticks(numpy.linspace(0, numpy.amax(data) * 1.1, 5))
                ax.set_ylabel(r'Production [$\text{ppm} / \text{s}$]')
                if scale == "long":
                    ax.set_xlabel("Time [days]")
                else:
                    ax.set_xlabel("Time [s]")
                plt.tight_layout()
                plt.savefig(f"figures/{path}-{k}-{eq:.1f}.pdf")
                plt.close()
    else:
        if isinstance(species, str):
            species = [species]
        if fig is None or ax is None:
            fig, ax = plt.subplots()
        fig.set_size_inches(6, 6.5)
        if entr_on:
            plt.axvspan(0, numpy.amax(time[time<=10]), facecolor=COLORS[-3], alpha=0.2, label="entrainment")
        for j, sp in enumerate(species):
            data = numpy.array(short_data[f"Y_{sp}"]) * 1e6
            if not mf:
                data *= mdot_exhaust / numpy.array(mass)
            mask = numpy.abs(data) <= limit
            for i in range(len(data)):
                data[i] = data[i] if not mask[i] else 0
            if normalized:
                data = data / numpy.amax(data)
            lb = sp.lower() if sp != "C5H8" else "isoprene"
            if markers:
                pltfcn(time, data, color=COLORS[j], label=lb, linestyle=linestyle, marker=markers[j])
            else:
                pltfcn(time, data, color=COLORS[j+CID], label=lb, linewidth=3, linestyle=linestyle)

        if legend_on:
            ax.legend(ncols=ncols, bbox_to_anchor=(0.5, 1.2), loc='upper center', fontsize=14)
        plt.xlabel("Time [s]")
        plt.ylabel(r'Production [$\text{ppm} / \text{s}$]')
        if mf:
            plt.yticks(numpy.logspace(numpy.log10(1e4), numpy.log10(1e6), 2))
        plt.grid(True)
        plt.tight_layout()
        if return_fig:
            return fig, ax
        # save after return
        if name is not None:
            plt.savefig(f"figures/{path}-{name}-{scale}-{eq:.1f}-{fp:.2f}.pdf", bbox_inches='tight')
        else:
            plt.savefig(f"figures/{path}-states-{scale}-{eq:.1f}-{fp:.2f}.pdf", bbox_inches='tight')

        plt.close()


def plot_total_HC_output_contour(path="minimal", name="", mode="hdf5", scale="long-term", model="atmosphere.yaml", mname="atmosphere", ekey="mdot_exhaust", rate=True):
    equiv_ratios = numpy.arange(0.2, 1.1, 0.1)
    farnesane_percents = numpy.arange(0, 0.21, 0.01)
    unformatted = f"{scale}-states-"+ "{:.1f}-{:.2f}" + f".{mode}"
    masses = numpy.zeros((len(farnesane_percents), len(equiv_ratios)))
    # get cantera solution
    atms_model = os.path.join(DATA_DIR, model)
    # creation of thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        atms = PlumeSolution(atms_model, name=mname, transport_model=None)
    # get all hydrocarbon species
    HC_species = []
    for sp in atms.species():
        if "H" in sp.composition.keys() and "C" in sp.composition.keys():
            HC_species.append(sp.name)
    for i, eq in enumerate(equiv_ratios):
        for j, fs in enumerate(farnesane_percents):
            if os.path.exists(os.path.join(path, unformatted.format(eq, fs))):
                # load yaml file for rate
                with open(os.path.join(path, f"thermo-states-{eq:.1f}-{fs:.2f}.yaml"), "r") as tsf:
                    yaml_data = yaml.load(tsf, Loader=yaml.SafeLoader)
                # load data
                if mode == "csv":
                    short_data = pd.read_csv(os.path.join(path, unformatted.format(eq, fs)), delimiter=",", engine="python")
                    for hc in HC_species:
                        masses[j][i] += short_data[f'Y_{hc}'].iloc[0] * yaml_data[ekey]
                elif mode == "hdf5":
                    short_data = h5py.File(os.path.join(path, unformatted.format(eq, fs)), "r")
                    for hc in HC_species:
                        if rate:
                            masses[j][i] += short_data[f'Y_{hc}'][0] * yaml_data[ekey] * 1000
                        else:
                            masses[j][i] += short_data[f'Y_{hc}'][0] * 1e6
            else:
                masses[j][i] = numpy.NAN
    # apply mask to masses
    mask = numpy.abs(masses) <= MASK_VALUE
    mask = numpy.logical_or(numpy.isnan(masses), mask)
    for i in range(len(equiv_ratios)):
        for j in range(len(farnesane_percents)):
            masses[j][i] = masses[j][i] if not mask[j][i] else MASK_VALUE/10
    # make contour
    X, Y = numpy.meshgrid(equiv_ratios, farnesane_percents)
    plt.contourf(X, Y, masses, cmap='viridis', norm=LogNorm(vmin=MASK_VALUE/10, vmax=numpy.amax(masses)))
    plt.xticks(numpy.linspace(0.2, 1.0, 5))
    plt.yticks(numpy.linspace(0, 0.2, 5))
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('Farnesane Fraction')
    cbar = plt.colorbar(label=r'Total Hydrocarbon Production [$\text{g} / \text{s}$]')
    cbar.ax.tick_params(labelsize=20)
    ticks = cbar.get_ticks()
    ticks = numpy.logspace(numpy.log10(ticks[0]), numpy.log10(numpy.amax(masses)), 6)
    tick_labels = ["$10^{" + f"{int(numpy.log10(i)):d}" + "}$" for i in ticks]
    # Modify only the first tick label
    tick_labels[0] = '$0$'
    # Set the modified tick labels
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(tick_labels)
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/{path}-hc-contour.pdf", bbox_inches='tight')
    plt.close()


def nox_comparison_plots(species=[], mode="hdf5", scale="long", eq=1.0, fp=0.10, name="", limit=None, pltfcn=plt.plot, ncols=2, abox=(0.5, 1.25)):
    # load data
    rates = []
    times = []
    hdata = []
    for p in ["ctp", "noxeq"]:
        with open(os.path.join(p, f"thermo-states-{eq:.1f}-{fp:.2f}.yaml"), "r") as tsf:
            mdot_exhaust = yaml.load(tsf, Loader=yaml.SafeLoader)["mdot_exhaust"]
        if mode == "hdf5":
            short_data = h5py.File(os.path.join(p, f'{scale}-term-states-{eq:.1f}-{fp:.2f}.hdf5'), "r")
            hdata.append(short_data)
            rates.append(mdot_exhaust / numpy.array(short_data["mass"]))
            times.append(numpy.array(short_data["time"]))

    if isinstance(species, str):
        species = [species]
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6.75)
    plt.axvspan(0, numpy.amax(times[0][times[0]<=10]), facecolor=COLORS[-3], alpha=0.2, label="entrainment")
    # Create some custom handles
    plt.plot([], [], color='k', linewidth=3, label="model")
    plt.plot([], [], color='k', linewidth=2, linestyle=":", label='equilibrium')
    for j, sp in enumerate(species):
        data0 = numpy.array(hdata[0][f"Y_{sp}"]) * rates[0] * 1e6
        data1 = numpy.array(hdata[1][f"Y_{sp}"]) * rates[1] * 1e6
        mask0 = numpy.abs(data0) <= MASK_VALUE
        mask1 = numpy.abs(data1) <= MASK_VALUE
        for i in range(len(data0)):
            data0[i] = data0[i] if not mask0[i] else 0
        for i in range(len(data1)):
            data1[i] = data1[i] if not mask1[i] else 0
        lb = sp.lower() if sp != "C5H8" else "isoprene"
        pltfcn(times[0], data0, color=COLORS[j], label=lb, linewidth=3)
        pltfcn(times[1], data1, color=COLORS[j], linestyle=":", linewidth=2)
    # handles, labels = plt.gca().get_legend_handles_labels()
    # handles += [handle1, handle2]
    plt.legend(ncols=ncols, bbox_to_anchor=abox, loc='upper center', fontsize=14)
    plt.xlabel("Time [s]")
    plt.ylabel(r'Production [$\text{ppm} / \text{s}$]')
    plt.grid(True)
    plt.tight_layout()
    if name:
        name = "-" + name
    plt.savefig(f"figures/nox-eq-{eq:.1f}-{fp:.2f}.pdf", bbox_inches='tight')
    plt.close()

def round_to_power_of_10(value):
    power = math.floor(math.log10(abs(value)))
    lower_power = 10 ** power
    upper_power = 10 ** (power + 1)
    return lower_power, upper_power

def make_water_states(sp, eq=0.9, path="water_ctp", mode="hdf5", scale="long", limit=MASK_VALUE):
    water_percents = numpy.arange(0, 0.05, 0.01)
    unformatted = f"{scale}-term-states-"+ "{:.1f}-0.10-{:.2f}" + f".{mode}"
    data = []
    times = []
    kept = []
    for i, wf in enumerate(water_percents):
        if os.path.exists(os.path.join(path, f"thermo-states-{eq:.1f}-0.10-{wf:.2f}.yaml")):
            kept.append(wf)
            with open(os.path.join(path, f"thermo-states-{eq:.1f}-0.10-{wf:.2f}.yaml"), "r") as tsf:
                mdot_exhaust = yaml.load(tsf, Loader=yaml.SafeLoader)["mdot_exhaust"]
            if mode == "hdf5":
                short_data = h5py.File(os.path.join(path, unformatted.format(eq, wf)), "r")
                times.append(numpy.array(short_data["time"]))
                data.append(numpy.array(short_data[f"Y_{sp}"]) * mdot_exhaust / numpy.array(short_data[f'mass']) * 1e6)
    # mask the data
    masked_data = []
    x = times[0]
    for d in range(len(data)):
        mask = numpy.abs(data[d]) <= limit
        temp_data = []
        for i in range(len(data[d])):
            if not mask[i]:
                temp_data.append(data[d][i])
            else:
                temp_data.append(limit / 10)
            # adjust data size
        ifc = interp1d(times[d], temp_data)
        masked_data.append(numpy.array(ifc(x)))
    masked_data = numpy.array(masked_data)
    # make plot
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6.5)
    plt.axvspan(0, numpy.amax(x[x<=10]), facecolor=COLORS[-3], alpha=0.2, label="entrainment")
    for i, r in enumerate(masked_data):
        plt.loglog(x, r, color=COLORS[i], label=f"{kept[i]:.2f}", linewidth=3)
    leg = plt.legend(bbox_to_anchor=(0.5, 1.32), loc='upper center', fontsize=16, ncols=3)
    leg.set_title("$X_{H2O}$")
    plt.xlabel('Time [s]')
    plt.ylabel(r'Production [$\text{ppm} / \text{s}$]')
    lb, ___ = round_to_power_of_10(numpy.amin(masked_data))
    ___, ub = round_to_power_of_10(numpy.amax(masked_data))
    if not os.path.exists("figures"):
        os.mkdir("figures")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/{path}-{sp.lower()}-{scale}-{eq:1.1f}.pdf", bbox_inches='tight')
    plt.close()



def entrainment_perturb_plot(species, path="entall", mode="hdf5", scale="long", eq=0.6, fp=0.10, pltfcn=plt.loglog):
    # load data
    with open(os.path.join(path, f"thermo-states-{eq:.1f}-{fp:.2f}-{1.0:.1f}.yaml"), "r") as tsf:
        mdot_exhaust = yaml.load(tsf, Loader=yaml.SafeLoader)["mdot_exhaust"]
    masses = []
    times = []
    entfs = [0.1, 1.0, 2.0, 5.0, 10.0]
    for v in entfs:
        if mode == "hdf5":
            short_data = h5py.File(os.path.join(path, f'{scale}-term-states-{eq:.1f}-{fp:.2f}-{v:.1f}.hdf5'), "r")
            masses.append(numpy.array(short_data[f"Y_{species}"]) * mdot_exhaust * 1e9 / numpy.array(short_data["mass"]))
            times.append(numpy.array(short_data["time"]))
    fig, ax = plt.subplots()
    plt.axvspan(0, numpy.amax(times[1][times[1]<=10]), facecolor=COLORS[-3], alpha=0.2, label="entrainment")
    for j, m in enumerate(masses):
        mask = numpy.abs(m) <= MASK_VALUE
        for i in range(len(m)):
            m[i] = m[i] if not mask[i] else 0
        pltfcn(times[j], m, color=COLORS[j], label=entfs[j], linewidth=3)
    leg = ax.legend()
    leg.set_title(r"$\dot{\omega}_{\text{ent}}$ scalar")
    plt.xlabel("Time [s]")
    plt.ylabel(r'Production [$\text{ppb} / \text{s}$]')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/{path}-entraincomp.pdf", bbox_inches='tight')
    plt.close()


def paper_plots():
    mode = "hdf5"
    path = "ctp"
    mapfile = os.path.join(DATA_DIR, "combustor_to_atmosphere.yaml")
    with open(mapfile, "r") as f:
        mapping = yaml.load(f, Loader=yaml.SafeLoader)
        mapping = {v: k for k, v in mapping.items()}
    all_specs = ["BENZENE", "TOLUENE",  "C5H8"]
    make_all = False
    eq1 = 0.9
    eq2 = 0.6
    for path in ["ctp", "ctpnn"]:
        if make_all:
                states_plots(None, path=path, mode=mode, scale="long", eq=eq1, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2)
                states_plots(None, path=path, mode=mode, scale="long", eq=eq2, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=3)
        else:
            states_plots(all_specs, path=path, mode=mode, scale="long", eq=eq1, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2)
            states_plots(all_specs, path=path, mode=mode, scale="long", eq=eq2, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2)
            states_plots(all_specs, path="rmin", mode=mode, scale="long", eq=eq1, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2)
            states_plots(all_specs, path="rmin", mode=mode, scale="long", eq=eq2, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2)
            nox_comparison_plots(all_specs, eq=eq1, fp=0.1, pltfcn=plt.loglog)
            nox_comparison_plots(all_specs, eq=eq2, fp=0.1, pltfcn=plt.loglog)
            states_plots(["O2", "N2"], path=path, mode=mode, scale="long", eq=eq2, fp=0.1, normalized=False, pltfcn=plt.loglog, name="o2", mf=True, ncols=3)
            states_plots(["O2", "N2"], path=path, mode=mode, scale="long", eq=eq1, fp=0.1, normalized=False, pltfcn=plt.loglog, name="o2", limit=1e-10, mf=True, ncols=3)
            for name in all_specs: #
                make_species_contour(name, index=0, path=path, mode=mode)
                make_species_contour(mapping[name], index=0, path="initials", mode=mode, scale="initial", ekey="mdot_in")
                make_water_states(name, path="water_ctp", eq=eq1, limit=1e-10)
                make_water_states(name, path="water_ctp", eq=eq2, limit=1e-10)
            make_temperature_contour(index=3, path=path, mode="yaml")
            make_temperature_contour(index=2, path=path, mode="yaml")
            # radical plots
            peroxys = ["BZBIPERO2", "TLBIPERO2", "ISOP34O2"] #"C6H5CO3"
            states_plots(peroxys, path=path, mode=mode, scale="long", eq=eq1, fp=0.1, name=f"peroxys-states", pltfcn=plt.loglog)
            states_plots(peroxys, path=path, mode=mode, scale="long", eq=eq2, fp=0.1, name=f"peroxys-states", pltfcn=plt.loglog)
            # ovocs = ["A2PANOO", "PHCOOH", "BZEMUCCO2H"] #"C6H5CO3"
            ovocs = ["BZBIPERNO3", "TLBIPERNO3",  "ISOP34NO3" ]
            states_plots(ovocs, path=path, mode=mode, scale="long", eq=eq1, fp=0.1, name=f"ovocs-states", pltfcn=plt.loglog)
            states_plots(ovocs, path=path, mode=mode, scale="long", eq=eq2, fp=0.1, name=f"ovocs-states", pltfcn=plt.loglog)
            # overall hydrocarbons
            HC_open = openAP_estimate()
            print(f"HC: {HC_open}")
            plot_total_HC_output_contour(path=path)
            plot_total_HC_output_contour(path="initials", scale="initial", mname="combustor", model="combustor.yaml", ekey="mdot_in")
            # entrainment perturb
            entrainment_perturb_plot("BENZENE")
            # surrogate speciation plots
            make_surrogate_bar_plots(all_specs, eq=eq1)
            make_surrogate_bar_plots(all_specs, eq=eq2)

def openAP_estimate():
    # A320 at cruise
    fuelflow = FuelFlow(ac='A320', eng='CFM56-5B4')
    emission = Emission(ac='A320', eng='CFM56-5B4')
    TAS = 493
    ALT = 35000
    FF = fuelflow.enroute(mass=48000, tas=TAS, alt=ALT)   # kg/s
    CO2 = emission.co2(FF)                    # g/s
    H2O = emission.h2o(FF)                    # g/s
    NOx = emission.nox(FF, tas=TAS, alt=ALT)  # g/s
    CO = emission.co(FF, tas=TAS, alt=ALT)    # g/s
    HC1 = emission.hc(FF, tas=TAS, alt=ALT)    # g/s

    # 737 idle
    fuelflow = FuelFlow(ac='B737', eng='CFM56-7B27')
    emission = Emission(ac='B737', eng='CFM56-7B27')
    TAS = 0
    ALT = 0
    FF = fuelflow.enroute(mass=48000, tas=TAS, alt=ALT)   # kg/s
    CO2 = emission.co2(FF)                    # g/s
    H2O = emission.h2o(FF)                    # g/s
    NOx = emission.nox(FF, tas=TAS, alt=ALT)  # g/s
    CO = emission.co(FF, tas=TAS, alt=ALT)    # g/s
    HC2 = emission.hc(FF, tas=TAS, alt=ALT)    # g/s
    return HC1, HC2


def update_phi_across():
    path = "initials"
    # create path to models
    fuel_model = os.path.join(DATA_DIR, "combustor.yaml")
    # creation of fuel thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, name="combustor", transport_model=None)
    for eq in numpy.arange(0.2, 1.1, 0.1):
        for fp in numpy.arange(0, 0.21, 0.01):
            # create initial fuel setting
            X_fuel = {"N-C10H22":0.4267, "I-C8H18":0.3302, "C7H8":0.2431}
            # adjust fuels for farnesane blend
            X_fuel = {k:v * (1-fp) for k,v in X_fuel.items()}
            # add farnesane blend
            Xfarne = 1 - numpy.sum([v for k,v in X_fuel.items()])
            if Xfarne > 0:
                X_fuel.update({"iC15H32": Xfarne})
            X_fuel = ", ".join([f"{k}:{v:.6f}" for k, v in X_fuel.items()])
            X_air = f"H2O: 0.01, O2:0.2095, N2:0.7808"
            multizone_combustor(fuel, 0.554, eq, X_fuel, X_air, norun=True, name_id=f"{eq:.1f}-{fp:.2f}", use_total_equiv=True)
            input()


def make_surrogate_bar_plots(species, path="surrogate", mode="hdf5", scale="long", eq=0.6, fp=0.10, pltfcn=plt.loglog):
    if isinstance(species, str):
        species = [species]
    all_data = []
    surrogates = [0, 2, 3, 5]
    for sn in surrogates:
        with open(os.path.join(path, f"thermo-states-{eq:.1f}-{fp:.2f}-{sn}.yaml"), "r") as tsf:
                mdot_exhaust = yaml.load(tsf, Loader=yaml.SafeLoader)["mdot_exhaust"]
        if mode == "hdf5":
            short_data = h5py.File(os.path.join(path, f'{scale}-term-states-{eq:.1f}-{fp:.2f}-{sn}.hdf5'), "r")
            mass = numpy.array(short_data["mass"])
        subdata = []
        for sp in species:
            value = numpy.array(short_data[f"Y_{sp}"][0]) * 1e6 * mdot_exhaust / numpy.array(mass)[0]
            subdata.append(value)
        all_data.append(subdata)
    all_data = numpy.array(all_data)
    all_data = numpy.transpose(all_data)
    # create bar plot
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6.5)
    bar_width = 0.35
    r1 = [i-bar_width for i in range(0, len(surrogates)*2, 2)]
    for j in range(len(species)):
        # Set the width of the bars
        values = [x + bar_width * j for x in r1]
        label = species[j].lower() if species[j].lower() != 'c5h8' else 'isoprene'
        ax.bar(values, all_data[j], color=COLORS[j], width=bar_width, edgecolor='k', label=label, hatch=HATCHES[j])
    # format figure
    ax.legend(ncols=3, bbox_to_anchor=(0.5, 1.1), loc='upper center', fontsize=16)
    ax.set_xticks([i for i in range(0, len(surrogates)*2, 2)])
    ax.set_xticklabels([f"$S_{i}$" for i in range(len(surrogates))])
    ax.set_ylabel(r'Production [$\text{ppm} / \text{s}$]', fontsize=20)
    ax.set_yscale("log")
    ax.yaxis.grid(True, color="k")
    plt.tight_layout()
    plt.savefig(f"figures/{path}-bar-{eq:.1f}-{fp:.2f}.pdf", bbox_inches='tight')
    plt.close()


def defense_plots():
    benz_sp = ["BENZENE", "BZBIPERO2", "BZBIPERNO3"]
    tol_sp = ["TOLUENE", "TLBIPERO2", "TLBIPERNO3"]
    iso_sp = ["C5H8", "ISOP34O2", "ISOP34NO3"]
    all_specs = ["BENZENE", "TOLUENE", "C5H8"]
    mode = "hdf5"
    for eq in [0.6, 0.9]:
        for path in ["ctp", "ctpnn"]:
            for n, sl in [("benz", benz_sp), ("tol", tol_sp), ("iso", iso_sp)]:
                states_plots(sl, path=path, mode=mode, scale="long", eq=eq, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2, name=n, legend_on=False)
                fig, ax = states_plots(sl, path=path, mode=mode, scale="long", eq=eq, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2, name=n, legend_on=False, return_fig=True, entr_on=False)
                states_plots(sl, path=path, mode=mode, scale="long", eq=eq, fp=0, normalized=False, pltfcn=plt.loglog, ncols=2, name=n, legend_on=False, fig=fig, ax=ax, linestyle=":", entr_on=False)
        for sp in all_specs:
            fig, ax = states_plots(sp, path="ctpnn", mode=mode, scale="long", eq=eq, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2, name=n, legend_on=False, return_fig=True, entr_on=True)
            fig, ax = states_plots(sp, path="ctp", mode=mode, scale="long", eq=eq, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2, name=n, legend_on=False, return_fig=True, entr_on=False, fig=fig, ax=ax, linestyle="-.", CID=1)
            fig, ax = states_plots(sp, path="ctp", mode=mode, scale="long", eq=eq, fp=0.0, normalized=False, pltfcn=plt.loglog, ncols=2, name=sp.lower(), legend_on=False, return_fig=True, entr_on=False, fig=fig, ax=ax, linestyle="--", CID=3)
            states_plots(sp, path="rmin", mode=mode, scale="long", eq=eq, fp=0.1, normalized=False, pltfcn=plt.loglog, ncols=2, name=sp.lower(), legend_on=False, entr_on=False, fig=fig, ax=ax, linestyle=":", CID=2)


if __name__ == "__main__":
    mode = "hdf5"
    path = "minimal"
    all_specs = ["TOLUENE", "BENZENE", "C5H8"]
    # paper_plots()
    defense_plots()
    # update_phi_across()
