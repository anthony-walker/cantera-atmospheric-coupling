import os
import re
import sys
import csv
import yaml
import h5py
import numpy
import random
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
    outdir = os.path.dirname(__file__)
    equiv_ratio, farnesane, sulfur = vals
    print(f"Running {equiv_ratio}, {farnesane}, {sulfur}...")
    log_name = f"combustor-{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}.log"
    log_name = os.path.join("combustor", log_name)
    log_file = open(log_name, "w")
    ostd = sys.stdout
    sys.stdout = log_file
    lf_closed = False
    try:
        combustor_atm_sim(equiv_ratio=equiv_ratio, nsteps=500, farnesane=farnesane, sulfur=sulfur, outdir=outdir)
    except Exception as e:
        lf_closed = True
        log_file.close()
        sys.stdout = ostd
        print(f"FAILED CASE: {equiv_ratio}, {farnesane}, {sulfur}...")
        case_files = filter(lambda x: f"{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}" in x, os.listdir("combustor"))
        for f in case_files:
            os.remove(os.path.join("combustor", f))
    # close if needed
    if not lf_closed:
        log_file.close()
        sys.stdout = ostd


def run_combustion_study(run_missing=True):
    equiv_ratios = numpy.arange(0.5, 2.75, 0.25)
    farnesane_percents = numpy.arange(0, 0.11, 0.02)
    sulfur_percents = [0.0001,]
    # get all combinations
    vals_list = []
    for eq in equiv_ratios:
        for xf in farnesane_percents:
            for xs in sulfur_percents:
                vals_list.append((eq, xf, xs))
    # filter to only missing
    if run_missing:
        all_files = filter(lambda x: "thermo-states" in x, os.listdir("combustor"))
        exists = [re.findall(r"\d+[.]\d+", f) for f in all_files]
        exists.sort()
        keep = []
        for eq, fs, sf in vals_list:
            case = [f"{eq:.1f}", f"{fs:.2f}", f"{sf:.4f}"]
            if case not in exists:
                keep.append((eq, fs, sf))
        vals_list = keep
    # random.shuffle(vals_list) # random shuffle of order
    # run multi-processing pool
    # equiv_ratio, nsteps, farnesane, sulfur, outdir
    with Pool(numpy.amin([os.cpu_count() - 2, len(vals_list)])) as p:
        p.map(wrapped, vals_list)


def convert_csv_to_hdf5(fname):
    print(f"Converting {fname} to hdf5...")
    with open(fname, "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        initial_values = [float(v) for v in next(reader)]
        for row in reader:
            vals = [float(v) for v in row]
            initial_values = numpy.vstack([initial_values, vals])
    # write hdf5 file
    hdf5_name = f"{fname.rsplit('.', maxsplit=1)[0]}.hdf5"
    hdf5_file = h5py.File(hdf5_name, "w")
    for i, h in enumerate(headers):
        hdf5_file.create_dataset(h, data=initial_values[:, i], dtype=numpy.double)
    hdf5_file.close()


def convert_all_files_to_hdf5():
    all_files = filter(lambda x: ".csv" in x and "atms" in x, os.listdir("combustor"))
    all_files = [os.path.join("combustor", f) for f in all_files]
    all_files = list(filter(lambda x: not os.path.isfile(f"{x.split('.')[0]}.hdf5"), all_files))
    with Pool(numpy.amin((os.cpu_count(), len(all_files)))) as p:
        p.map(convert_csv_to_hdf5, all_files)


def get_hdf5_data(hf_name, xkey, ykey):
    with h5py.File(hf_name, "r") as hf:
        x = hf[xkey][:]
        y = hf[ykey][:]
    return x, y

def get_yaml_data(k1, k2, oi1=None, oi2=None, ffcn=None, skey=None,):
    af = filter(lambda x: "thermo-states" in x, os.listdir("combustor/yamls"))
    x = []
    y = []
    eq_ratio = []
    farn = []
    sulf = []
    for f in af:
        with open(os.path.join("combustor/yamls", f), "r") as cf:
            data = yaml.safe_load(cf)
        cx = data[k1]
        if oi1 is not None:
            cx = cx[oi1]
        cy = data[k2]
        if oi2 is not None:
            cy = cy[oi2]
        x.append(cx)
        y.append(cy)
        eq_ratio.append(data["equivalence_ratio"])
        farn.append(data["farnesane"])
        sulf.append(data["sulfur"])
    data = list(zip(x, y, eq_ratio, farn, sulf))
    if ffcn is not None:
        data = list(filter(ffcn, data))
    x, y, eq_ratio, farn, sulf = zip(*data)
    x = [float(q) for q in x]
    y = [float(q) for q in y]
    if skey is not None:
        xy_data = list(zip(x, y))
        xy_data.sort(key=skey)
        x, y = zip(*xy_data)
    # normalize
    return x, y

def plot_temperature_equiv_ratio():
    # x, y, eq_ratio, farn, sulf - data order
    x1, y1 = get_yaml_data("equivalence_ratio", "T", oi2=3, ffcn=lambda x: x[3] == "0.10", skey=lambda x: x[0])
    plt.plot(x1, y1)
    plt.show()
    plt.close()

def plot_mass_time(speckey, name, ylab):
    # figure path
    fpath = f"combustor/figures/{name}"
    if not os.path.isdir(fpath):
        os.mkdir(fpath)
    # get file name by arguments
    ghfn = lambda eq, fs, sf: os.path.join("combustor", f"atms-states-{eq:.1f}-{fs:.2f}-{sf:.4f}.hdf5")
    gethf5 = lambda x: os.path.join("combustor", x)
    # get all files
    all_files = os.listdir("combustor")
    # all cases
    equiv_ratios = numpy.arange(0.5, 2.75, 0.25)
    farnesane_percents = numpy.arange(0, 0.11, 0.02)
    sulfur_percents = [0.0001,]
    # isoprene
    for phi in equiv_ratios:
        rfiles = filter(lambda x: f"atms-states-{phi:.1f}" in x and x.endswith(".hdf5"), all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), speckey, "time")
            nc, tc = get_hdf5_data(gethf5(rf), "mass", "time")
            x_iso.append(tc)
            y_iso.append(xc * nc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Time [s]", fontsize=18)
        ax.set_ylabel(ylab, fontsize=22)
        ax.set_xlim([0, 0.05])
        y0 = 10 ** numpy.floor(numpy.log10(y_min)) if y_min > 0 else y_min
        yf = 10 ** numpy.ceil(numpy.log10(y_max)) if y_max > 0 else y_max
        # ax.set_ylim([y0, yf])
        ax.set_xticks(numpy.arange(0, 0.06, 0.01))
        llab = "$X_{{iC15H32}}: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[1]
            ax.plot(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.45, 1.03))
        fig.savefig(os.path.join(fpath, f"{name}-mass-{phi:.1f}.pdf"), bbox_inches='tight')
        plt.close()

    for phi in farnesane_percents:
        rfiles = filter(lambda x: f"-{phi:.2f}-" in x and x.endswith(".hdf5") and "atms" in x, all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), speckey, "time")
            nc, tc = get_hdf5_data(gethf5(rf), "mass", "time")
            x_iso.append(tc)
            y_iso.append(xc * nc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Time [s]", fontsize=18)
        ax.set_ylabel(ylab, fontsize=22)
        ax.set_xlim([0, 0.05])
        y0 = 10 ** numpy.floor(numpy.log10(y_min)) if y_min > 0 else y_min
        yf = 10 ** numpy.ceil(numpy.log10(y_max)) if y_max > 0 else y_max
        # ax.set_ylim([y0, yf])
        ax.set_xticks(numpy.arange(0, 0.06, 0.01))
        llab = "$\phi: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[0]
            ax.plot(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.3, 1.03))
        fig.savefig(os.path.join(fpath, f"{name}-mass-{phi:.2f}.pdf"), bbox_inches='tight')
        plt.close()

def plot_species_time(speckey, name, ylab, auto=False):
    # figure path
    fpath = f"combustor/figures/{name}"
    if not os.path.isdir(fpath):
        os.mkdir(fpath)
    # get file name by arguments
    ghfn = lambda eq, fs, sf: os.path.join("combustor", f"atms-states-{eq:.1f}-{fs:.2f}-{sf:.4f}.hdf5")
    gethf5 = lambda x: os.path.join("combustor", x)
    # get all files
    all_files = os.listdir("combustor")
    # all cases
    equiv_ratios = numpy.arange(0.5, 2.75, 0.25)
    farnesane_percents = numpy.arange(0, 0.11, 0.02)
    sulfur_percents = [0.0001,]
    # isoprene
    for phi in equiv_ratios:
        rfiles = filter(lambda x: f"atms-states-{phi:.1f}" in x and x.endswith(".hdf5"), all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), speckey, "time")
            x_iso.append(tc)
            y_iso.append(xc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Time [s]", fontsize=18)
        ax.set_ylabel(ylab, fontsize=22)
        ax.set_xlim([0, 0.05])
        if not auto:
            y0 = 10 ** numpy.floor(numpy.log10(y_min)) if y_min > 0 else y_min
            yf = 10 ** numpy.ceil(numpy.log10(y_max)) if y_max > 0 else y_max
            ax.set_ylim([y0, yf])
        ax.set_xticks(numpy.arange(0, 0.06, 0.01))
        llab = "$X_{{iC15H32}}: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[1]
            ax.semilogy(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.45, 1.03))
        fig.savefig(os.path.join(fpath, f"{name}-{phi:.1f}.pdf"), bbox_inches='tight')
        plt.close()

    for phi in farnesane_percents:
        rfiles = filter(lambda x: f"-{phi:.2f}-" in x and x.endswith(".hdf5") and "atms" in x, all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), speckey, "time")
            x_iso.append(tc)
            y_iso.append(xc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Time [s]", fontsize=18)
        ax.set_ylabel(ylab, fontsize=22)
        ax.set_xlim([0, 0.05])
        if not auto:
            y0 = 10 ** numpy.floor(numpy.log10(y_min)) if y_min > 0 else y_min
            yf = 10 ** numpy.ceil(numpy.log10(y_max)) if y_max > 0 else y_max
            ax.set_ylim([y0, yf])
        ax.set_xticks(numpy.arange(0, 0.06, 0.01))
        llab = "$\phi: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[0]
            ax.semilogy(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.3, 1.03))
        fig.savefig(os.path.join(fpath, f"{name}-{phi:.2f}.pdf"), bbox_inches='tight')
        plt.close()

def plot_temperature_time(auto=False):
    # figure path
    fpath = f"combustor/figures/temperature"
    if not os.path.isdir(fpath):
        os.mkdir(fpath)
    # get file name by arguments
    ghfn = lambda eq, fs, sf: os.path.join("combustor", f"atms-states-{eq:.1f}-{fs:.2f}-{sf:.4f}.hdf5")
    gethf5 = lambda x: os.path.join("combustor", x)
    # get all files
    all_files = os.listdir("combustor")
    # all cases
    equiv_ratios = numpy.arange(0.5, 2.75, 0.25)
    farnesane_percents = numpy.arange(0, 0.11, 0.02)
    sulfur_percents = [0.0001,]
    # isoprene
    for phi in equiv_ratios:
        rfiles = filter(lambda x: f"atms-states-{phi:.1f}" in x and x.endswith(".hdf5"), all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 0
        y_min = 10000
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), "T", "time")
            x_iso.append(tc)
            y_iso.append(xc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Time [s]", fontsize=18)
        ax.set_ylabel("Temperature [K]", fontsize=18)
        ax.set_xlim([0, 0.05])
        if not auto:
            y0 = int(y_min / 1.1)
            yf = int(y_max * 1.1)
            ax.set_ylim([y0, yf])
        ax.set_xticks(numpy.arange(0, 0.06, 0.01))
        llab = "$X_{{iC15H32}}: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[1]
            ax.plot(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.45, 1.03))
        fig.savefig(os.path.join(fpath, f"temperature-{phi:.1f}.pdf"), bbox_inches='tight')
        plt.close()

    for phi in farnesane_percents:
        rfiles = filter(lambda x: f"-{phi:.2f}-" in x and x.endswith(".hdf5") and "atms" in x, all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), "T", "time")
            x_iso.append(tc)
            y_iso.append(xc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Time [s]", fontsize=18)
        ax.set_ylabel("Temperature [K]", fontsize=18)
        ax.set_xlim([0, 0.05])
        if not auto:
            y0 = int(y_min / 1.1)
            yf = int(y_max * 1.1)
            ax.set_ylim([y0, yf])
        ax.set_xticks(numpy.arange(0, 0.06, 0.01))
        llab = "$\phi: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[0]
            ax.plot(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.3, 1.03))
        fig.savefig(os.path.join(fpath, f"temperature-{phi:.2f}.pdf"), bbox_inches='tight')
        plt.close()

def plot_species_temperature(speckey, name, ylab, auto=False):
    # figure path
    fpath = f"combustor/figures/{name}"
    if not os.path.isdir(fpath):
        os.mkdir(fpath)
    # get file name by arguments
    ghfn = lambda eq, fs, sf: os.path.join("combustor", f"atms-states-{eq:.1f}-{fs:.2f}-{sf:.4f}.hdf5")
    gethf5 = lambda x: os.path.join("combustor", x)
    # get all files
    all_files = os.listdir("combustor")
    # all cases
    equiv_ratios = numpy.arange(0.5, 2.75, 0.25)
    farnesane_percents = numpy.arange(0, 0.11, 0.02)
    sulfur_percents = [0.0001,]
    # isoprene
    for phi in equiv_ratios:
        rfiles = filter(lambda x: f"atms-states-{phi:.1f}" in x and x.endswith(".hdf5"), all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), speckey, "T")
            x_iso.append(tc)
            y_iso.append(xc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Temperature [K]", fontsize=18)
        ax.set_ylabel(ylab, fontsize=22)
        if not auto:
            y0 = 10 ** numpy.floor(numpy.log10(y_min)) if y_min > 0 else y_min
            yf = 10 ** numpy.ceil(numpy.log10(y_max)) if y_max > 0 else y_max
            ax.set_ylim([y0, yf])
        llab = "$X_{{iC15H32}}: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[1]
            ax.semilogy(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.45, 1.03))
        fig.savefig(os.path.join(fpath, f"{name}-{phi:.1f}-temp.pdf"), bbox_inches='tight')
        plt.close()

    for phi in farnesane_percents:
        rfiles = filter(lambda x: f"-{phi:.2f}-" in x and x.endswith(".hdf5") and "atms" in x, all_files)
        x_iso = []
        y_iso = []
        # sort files
        rfiles = list(rfiles)
        rfiles.sort()
        y_max = 1e-16
        y_min = 1
        for rf in rfiles:
            xc, tc = get_hdf5_data(gethf5(rf), speckey, "T")
            x_iso.append(tc)
            y_iso.append(xc)
            y_min = numpy.amin(xc) if numpy.amin(xc) < y_min else y_min
            y_max = numpy.amax(xc) if numpy.amax(xc) > y_max else y_max
        ct = 0
        # create plot and figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_xlabel("Temperature [K]", fontsize=18)
        ax.set_ylabel(ylab, fontsize=22)
        if not auto:
            y0 = 10 ** numpy.floor(numpy.log10(y_min)) if y_min > 0 else y_min
            yf = 10 ** numpy.ceil(numpy.log10(y_max)) if y_max > 0 else y_max
            ax.set_ylim([y0, yf])
        llab = "$\phi: {:s}$"
        for cx, cy in zip(x_iso, y_iso):
            fsn = re.findall(r'\d+.\d+', rfiles[ct])[0]
            ax.semilogy(cx, cy, color=cbcolors[ct], linewidth=2, label=llab.format(fsn))
            ct += 1
        ax.legend(ncols=1, bbox_to_anchor=(1.3, 1.03))
        fig.savefig(os.path.join(fpath, f"{name}-{phi:.2f}-temp.pdf"), bbox_inches='tight')
        plt.close()


def heat_release_plots():
    # Plot results
    f, ax1 = plt.subplots(1, 1)
    ax1.plot(states.tres, states.heat_release_rate, '.-', color='C0')
    ax2 = ax1.twinx()
    ax2.plot(states.tres[:-1], states.T[:-1], '.-', color='C1')
    ax1.set_xlabel('residence time [s]')
    ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
    ax2.set_ylabel('temperature [K]', color='C1')
    f.tight_layout()
    plt.savefig(os.path.join(outdir, "combustor", f"heat-release-rate-{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}.pdf"))
    plt.close()


def make_plots():
    configs = [("Y_C5H8", "isoprene", "$Y_{C5H8}$"), ("Y_C7H8", "toluene", "$Y_{C7H8}$"), ("Y_A1", "benzene", "$Y_{C6H6}$"), ("Y_O2", "oxygen", "$Y_{O2}$"), ("Y_O3", "ozone", "$Y_{O3}$"), ("Y_N2", "nitrogen", "$Y_{N2}$"), ("Y_NO2", "no2", "$Y_{NO2}$"), ("Y_NO3", "no3", "$Y_{NO3}$"), ("Y_NO", "no", "$Y_{NO}$"), ("Y_CO2", "co2", "$Y_{CO2}$"), ("Y_CO", "co", "$Y_{CO}$"), ("Y_OH", "oh", "$Y_{OH}$")]
    for x, n, xl in configs:
        try:
            plot_species_time(x, n, xl, auto=True)
            plot_species_temperature(x, n, xl, auto=True)
            # plot_mass_time(x, n, xl)
        except Exception as e:
            print(f"Failed for {x}...")
    plot_temperature_time(auto=True)

if __name__ == "__main__":
    # combustor_atm_sim(0.6, 500, 0.08, 0.0001, "./")
    # run_combustion_study()
    make_plots()
    # plot_species_time(*("Y_C5H8", "isoprene", "$Y_{C5H8}$"))
    # plot_species_time(*("Y_C5H8", "isoprene", "$Y_{C5H8}$"), auto=True)

