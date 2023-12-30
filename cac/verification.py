import os
import re
import csv
import h5py
import yaml
import numpy
import click
import warnings
import pandas as pd
import cantera as ct
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from cac.constants import DATA_DIR, COLORS
from cac.combustor import multizone_combustor
from cac.combustor import curve_fit_thrust_data, reformat_hdf5
from cac.reactors import PlumeSolution, PlumeReactor


def mixing_box_model_verification():
    # setup plume reactor
    plume_gas = PlumeSolution('gri30.yaml')
    T_amb = 227
    plume_gas.TP = T_amb, ct.one_atm
    plume_gas.set_equivalence_ratio(0, 'CH4:1.0', 'O2:1.0, N2:3.76')
    air_state = plume_gas.state # get the state of air
    air_enthalpy = plume_gas.partial_molar_enthalpies
    air_cp = plume_gas.partial_molar_cp
    plume_gas.TPX = 1000, 0.2 * ct.one_atm, "H2O:1.0, O2:0.5, N2:3.76"
    # plume_gas.TP = 400, 0.2 * ct.one_atm
    # create plume reactor
    pr = PlumeReactor(plume_gas)
    pr.volume = 1.0
    pr.TX_air = [air_state[0]] + list(air_state[2:])
    pr.enthalpy_air = air_enthalpy
    pr.cp_air = air_cp
    # create network
    net = ct.ReactorNet([pr])
    net.preconditioner = ct.AdaptivePreconditioner()
    pstates = ct.SolutionArray(plume_gas)
    times = []
    # create data to plot
    while net.time < 10e4:
        pstates.append(pr.thermo.state)
        times.append(net.time)
        for i in range(10):
            net.step()
    # final state
    pstates.append(pr.thermo.state)
    times.append(net.time)
    # plot verification plots
    T0 = 1000
    normalized_temp = [(T - T_amb) / (T0 - T_amb) for T in pstates.T]
    h2o_id = plume_gas.species_index("H2O")
    Xh2o = [arr[h2o_id] for arr in pstates.X]
    normalized_Xh2o = [i/Xh2o[0] for i in Xh2o]
    # plot data
    plt.loglog(times, normalized_temp, color=COLORS[0], linestyle="-.", linewidth=2, label="$T_{\\text{ct}}$")
    plt.loglog(times, normalized_Xh2o, color=COLORS[1], linestyle="--", linewidth=2, label="$X_{\\text{H2O, ct}}$")
    plt.xlim([10e-4, 10])
    plt.ylim([0.01, 10])
    # load data
    with open(os.path.join(DATA_DIR, "verification", "box-model-verification.csv"), "r") as f:
        rdr = csv.reader(f)
        next(rdr)
        next(rdr)
        data = zip(*[r for r in rdr])
        xT, yT, xX, yX = [list(map(float, filter(lambda x: x, l))) for l in data]
    plt.loglog(xT, yT, color=COLORS[0], linestyle="", marker="s", label="$T_{\\text{Karcher}}$")
    plt.loglog(xX, yX, color=COLORS[1], linestyle="", marker="s", label="$X_{\\text{H2O, Karcher}}$")
    plt.ylabel("Normalized variables")
    plt.xlabel("Time [s]")
    plt.legend()
    fname = os.path.join(DATA_DIR, "verification", "box-model-verification.pdf")
    # add RMS error file
    plt.savefig(fname)


def combustor_design_params():
    # design parameters
    equiv_ratio = 1.77
    EPR = 25.8 # engine pressure ratio
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    X_air = "O2:1.0, N2:3.76"
    X_fuel = ", ".join([f'{k}:{v}'for k, v in {"NC10H22":0.4267, "IC8H18":0.3302, "C7H8":0.2431}.items()])
    X_fuel = "POSF10325:1.0"
    # create path to fuel model
    fuel_model = os.path.join(DATA_DIR, "A2NOx.yaml")
    # fuel_model = "h2o2.yaml"
    # X_fuel = "H2:1.0"
    # creation of fuel thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, transport_model=None, basis="mole")
    fuel.TPX = T_atm, p_atm, X_air
    # isentropic compression
    p1 = p_atm * EPR
    T1 = T_atm * (p1 / p_atm) ** (ct.gas_constant / fuel.cp_mole)
    # set input fuel conditions
    fuel.TP = T1, p1
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, X_air, basis="mole")
    return fuel, equiv_ratio, X_fuel, X_air


def parallel_run_combustor(x):
    fuel, equiv_ratio, X_fuel, X_air = combustor_design_params()
    out_dir = ver_dir = os.path.join(DATA_DIR, "verification")
    # try:
    print(f"Running thrust-level: {x:.3f}")
    combustor, thermo_states = multizone_combustor(fuel, x, equiv_ratio, X_fuel, X_air, n_pz=21, name_id=f"-verification-{x:0.2f}", outdir=out_dir)
    # except Exception as e:
    #     print(f"Failed for Thrust-level: {x:.3f}")
    #     # print(e)
    #     return None
    return (combustor.thermo.state, combustor.mass, thermo_states, x)


@click.command()
@click.option('--regenerate', is_flag=True, help='Regenerate verification data')
def combustor_verification(regenerate):
    # curve_fit_thrust_data(test=True)
    ver_dir = os.path.join(DATA_DIR, "verification")
    cbs = os.path.join(ver_dir, f"combustor-verification.hdf5")
    yfile = os.path.join(ver_dir, "combustor-verification.yaml")
    if not os.path.isdir(ver_dir):
        os.mkdir(ver_dir)
    fuel, equiv_ratio, X_fuel, X_air = combustor_design_params()
    thrust_levels = numpy.linspace(0.1, 1.0, 10)
    thrust_levels = [1.0]
    # thrust_levels = numpy.array([0.07, 0.3, 0.65, 0.85, 1.0])
    if not os.path.isfile(cbs) or regenerate:
        # run in parallel
        with mp.Pool(os.cpu_count()) as pool:
            res = pool.map(parallel_run_combustor, thrust_levels)
        # solution array
        combustor_states = ct.SolutionArray(fuel, extra=["mass", "TL", "phi", "mdot", "fuel_mass"])
        verification_data = {}
        for r in res:
            if r is not None:
                state, mass, thermo_states, tl = r
                verification_data[f"thrust-level-{tl}"] = thermo_states
                combustor_states.append(state, mass=mass, TL=tl, phi=equiv_ratio, mdot=1.086*tl, fuel_mass=thermo_states["initial_fuel_mass"])
        # dump verification data
        os.path.join("")
        with open(yfile, "w") as f:
            yaml.dump(verification_data, f)
        # write out state and reformats
        combustor_states.save(cbs, overwrite=True, name="thermo")
    # create plots of profiles from yaml
    with open(yfile, "r") as f:
        yaml_data = yaml.load(f, Loader=yaml.SafeLoader)
    # temperature profiles
    fig, ax = plt.subplots()
    for i, k in enumerate(yaml_data.keys()):
        ldata = yaml_data[k]
        ax.plot(ldata["phi_distribution"], ldata["temperature_distribution"], color=COLORS[i], label=f"{float(k.split('-')[-1]):.2f}")
    ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_ylabel("Temperature [K]")
    ax.set_xlabel("Equivalence Ratio")
    fig.savefig(os.path.join(ver_dir, f"T-profiles.pdf"), bbox_inches='tight')
    plt.close()

    # mass flow fraction profiles
    fig, ax = plt.subplots()
    for i, k in enumerate(yaml_data.keys()):
        ldata = yaml_data[k]
        ax.plot(ldata["phi_distribution"], ldata["mass_flow_fraction_distribution"], color=COLORS[i], label=f"{float(k.split('-')[-1]):.2f}")
    ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_ylabel("Fraction of Total Mass Flow")
    ax.set_xlabel("Equivalence Ratio")
    fig.savefig(os.path.join(ver_dir, f"mass-flow-fraction-profiles.pdf"), bbox_inches='tight')
    plt.close()

    # create plot with generated data
    with h5py.File(cbs, "r") as hf:
        Y_data = hf["thermo"]["data"]["Y"]
        mass_data = hf["thermo"]["data"]["mass"]
        fuel_mass_data = hf["thermo"]["data"]["fuel_mass"]
        mdot_data = hf["thermo"]["data"]["mdot"]
        thrust_levels = numpy.array(hf["thermo"]["data"]["TL"])
        # compute NOx
        mass_nox = Y_data[:, fuel.species_index("NO")] * mass_data * 1000 # grams
        nox_ratio = mass_nox / fuel_mass_data
        mass_co = Y_data[:, fuel.species_index("CO")] * mass_data * 1000 # grams
        co_ratio = mass_co / fuel_mass_data
    # plots
    fig, axes = plt.subplots(1, 2)
    with open(os.path.join(ver_dir, "NOX-data.csv")) as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)
        nox_x, nox_y = zip(*list(map(lambda x: (float(x[0])/100, float(x[1])), reader)))
    # plot
    axes[0].plot(thrust_levels[::-1], nox_ratio, color=COLORS[0], label="Model prediction")
    axes[0].plot(nox_x, nox_y, color=COLORS[1], linestyle="", marker="o", label="EDB data")
    axes[0].set_xlabel("Mass Flow Fraction")
    axes[0].set_ylabel("EI $\\text{NO}_x$ [g/$\\text{kg}_{\\text{fuel}}$]")

    with open(os.path.join(ver_dir, "CO-data.csv")) as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)
        co_x, co_y = zip(*list(map(lambda x: (float(x[0])/100, float(x[1])), reader)))
    # plot
    axes[1].plot(thrust_levels[::-1], co_ratio, color=COLORS[0], label="Model prediction")
    axes[1].plot(co_x, co_y, color=COLORS[1], linestyle="", marker="o", label="EDB data")
    axes[1].set_xlabel("Mass Flow Fraction")
    axes[1].set_ylabel("EI CO [g/$\\text{kg}_{\\text{fuel}}$]")
    plt.subplots_adjust(wspace=0.25)
    fig.savefig(os.path.join(ver_dir, f"emissions-indices.pdf"), bbox_inches='tight')
    plt.close()

    # temperature as a function of z
    # fast mixing data
    fmh5 = os.path.join(ver_dir, "fast-mixing-states-verification-1.00.hdf5")
    with h5py.File(fmh5, "r") as fhf:
        z_fast_data = numpy.array(fhf["z"])
        z_fast_data /= numpy.amax(z_fast_data)
        z_fast_data *= 100
        T_fast_data = numpy.array(fhf["T"])
    # slow mixing data
    smh5 = os.path.join(ver_dir, "slow-mixing-states-verification-1.00.hdf5")
    with h5py.File(smh5, "r") as shf:
        z_slow_data = numpy.array(shf["z"])
        z_slow_data /= numpy.amax(z_slow_data)
        z_slow_data *= 100
        T_slow_data = numpy.array(shf["T"])
    # fast mixing data
    fmh5 = os.path.join(ver_dir, "fast-mixing-states-verification-0.07.hdf5")
    with h5py.File(fmh5, "r") as fhf7:
        z_fast_data_s = numpy.array(fhf7["z"])
        z_fast_data_s /= numpy.amax(z_fast_data_s)
        z_fast_data_s *= 100
        T_fast_data_s = numpy.array(fhf7["T"])
    # slow mixing data
    smh5 = os.path.join(ver_dir, "slow-mixing-states-verification-0.07.hdf5")
    with h5py.File(smh5, "r") as shf7:
        z_slow_data_s = numpy.array(shf7["z"])
        z_slow_data_s /= numpy.amax(z_slow_data_s)
        z_slow_data_s *= 100
        T_slow_data_s = numpy.array(shf7["T"])
    # csv data
    with open (os.path.join(ver_dir, "sz-temperatures.csv")) as f:
        rdr = csv.reader(f)
        next(rdr)
        next(rdr)
        zfv, Tfv, zsv, Tsv = [list(map(float, filter(lambda x: x, k))) for k in zip(*[r for r in rdr ])]
    fig, ax = plt.subplots()
    ax.plot(z_fast_data, T_fast_data, color=COLORS[0], label="fast - 100%")
    ax.plot(z_slow_data, T_slow_data, color=COLORS[1], label="slow - 100%")
    ax.plot(z_fast_data_s, T_fast_data_s, color=COLORS[0], label="fast - 7%")
    ax.plot(z_slow_data_s, T_slow_data_s, color=COLORS[1], label="slow - 7 %")
    ax.plot(zfv, Tfv, color=COLORS[0], linestyle="", marker="o", label="fast source")
    ax.plot(zsv, Tsv, color=COLORS[1], linestyle="", marker="o", label="slow source")

    ax.set_xlim([0, 100.0])
    ax.set_xlabel("Length [%]")
    ax.set_ylabel("Temperature [K]")
    ax.legend()

    plt.show()

def convert_mission_out():
    ver_dir = os.path.join(DATA_DIR, "verification")
    fname = os.path.join(ver_dir, "CFM56_5B_737Mission.out")
    with open(fname, "r") as f:
        lines = f.read().split("\n")
        headers = ", ".join(re.split("\s+", lines[2])[1:-1])
        nls = [headers]
        for l in lines[3:]:
            r = re.split("\s+", l)[2:-1]
            nls.append(", ".join(r))
        content = "\n".join(nls)
    with open(os.path.join(ver_dir, "CFM56-5B-737.csv"), "w") as f:
        f.write(content)

def check_thrust():
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", 'CFM56-7B-737.csv'), delimiter=",", engine="python")
    ma = cycle_data['W3[kg/s]'].values
    mf = cycle_data['Wf[kg/s]'].values
    thrust = cycle_data['NetThrust[kN]'].values
    for i, j, k in zip(ma, mf, thrust):
        print(i, j, k)

def mcm_verification():
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    volume = 1e3
    gas.TPX = 298, ct.one_atm, "O3:0.2, CO:0.2, O2:0.2, N2:0.2, H2O:0.2"
    # setup initial conditions as 30 ppbv O3 and 100 ppbv CO with 1% water vapor
    r = PlumeReactor(gas)
    r.entrainment = False # Turn off entrainment
    r.volume = volume
    # setup reactor network
    net = ct.ReactorNet([r])
    net.max_steps = 1e7
    net.preconditioner = ct.AdaptivePreconditioner()
    # set tolerances
    # net.rtol = 1e-10
    # net.atol = 1e-16
    # setup time steps
    nsteps = 10
    ndays = 3
    seconds = ndays * 24 * 3600 # convert days to seconds
    steps = numpy.linspace(0, seconds, nsteps)
    # setup solution array
    arr = ct.SolutionArray(gas, extra=["t", "mass"])
    arr.append(r.thermo.state, t=0.0, mass=r.mass)
    # loop through all time steps
    for t in numpy.linspace(0, seconds, 500): #for t in steps:
        loop = True
        while loop:
            try:
                print(f"Advancing to {t:0.2f}..")
                net.advance(t)
                net.max_time_step = 0
                loop = False
            except Exception as e:
                print("Reducing max time step")
                if net.max_time_step == 0:
                    net.max_time_step == 1
                else:
                    net.max_time_step = net.max_time_step / 5
        # append the state
        arr.append(r.thermo.state, t=net.time, mass=r.mass)
    # save hdf5 file
    hdf5_file = os.path.join(ver_dir, "mcm.hdf5")
    arr.save(hdf5_file, overwrite=True, name="thermo")
    reformat_hdf5(hdf5_file)
    # verification data
    mcm_data = pd.read_csv(os.path.join(ver_dir, 'mcm-verification.csv'), engine="python")
    mcm_time = mcm_data['TIME'].values
    # print(gas.reactions()[0].rate_coeff_units)
    # plot from arr object
    fig, ax = plt.subplots()
    mcm_o3 = mcm_data["O3"].values
    o3_id = gas.species_index("O3")
    ax.plot(arr.t / 24 / 3600, arr.Y[:, o3_id] * arr.D / gas.molecular_weights[o3_id] * 1e9, color=COLORS[1])
    ax.set_xlim(0, 3)
    ax.set_xticks(numpy.linspace(0, 3, 4))
    # ax.set_xticklabels(labels)
    # ax.plot(mcm_time, mcm_o3, color=COLORS[0], linestyle="", marker="o")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("$O_3$")
    plt.show()



if __name__ == "__main__":
    # get_mass_flow_thrust_data(0.5)
    check_thrust()
