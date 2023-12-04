import os
import re
import sys
import csv
import h5py
import numpy
import random
import matplotlib.pyplot as plt
from cac.combustor import combustor_atm_sim
from multiprocessing import Pool

# change font of plots
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'


def wrapped(vals):
    outdir = os.path.dirname(__file__)
    equiv_ratio, farnesane, sulfur = vals
    print(f"Running {equiv_ratio}, {farnesane}, {sulfur}...")
    log_name = f"combustor-{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}.log"
    log_name = os.path.join("combustor", log_name)
    log_file = open(log_name, "w")
    ostd = sys.stdout
    sys.stdout = log_file
    combustor_atm_sim(equiv_ratio=equiv_ratio, nsteps=500, farnesane=farnesane, sulfur=sulfur, outdir=outdir)
    log_file.close()
    sys.stdout = ostd


def run_combustion_study(run_missing=True):
    equiv_ratios = numpy.arange(0.5, 2.6, 0.1)
    farnesane_percents = numpy.arange(0, 0.11, 0.01)
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
    random.shuffle(vals_list)
    # run multi-processing pool
    # equiv_ratio, nsteps, farnesane, sulfur, outdir
    with Pool(numpy.amin([os.cpu_count() - 4, len(vals_list)])) as p:
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
    all_files = filter(lambda x: ".csv" in x, os.listdir("combustor"))
    all_files = [os.path.join("combustor", f) for f in all_files]
    all_files = list(filter(lambda x: not os.path.isfile(f"{x.split('.')[0]}.hdf5"), all_files))
    with Pool(numpy.amin((os.cpu_count() - 2, len(all_files)))) as p:
        p.map(convert_csv_to_hdf5, all_files)


def plot_temperature(hf_name):
    with h5py.File(hf_name, "r") as hf:
        temp = hf['X_CH*'][:]
        time = hf['time'][:]
    plt.plot(time, temp)
    plt.show()


if __name__ == "__main__":
    # combustor_atm_sim(0.6, 500, 0.08, 0.0001, "./")
    run_combustion_study()
    convert_all_files_to_hdf5()
    # plot_temperature("combustor/atms-states-0.7-0.05-0.0001.hdf5")
