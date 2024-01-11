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
from scipy.optimize import fsolve
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from cac.constants import DATA_DIR, COLORS
from cac.combustor import multizone_combustor, get_prop_by_thrust_level
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
    # X_fuel = "NC12H26:0.404, IC8H18:0.295, TMBENZ:0.073, NPBENZ:0.228"
    # X_fuel = ", ".join([f'{k}:{v}'for k, v in {"N-C10H22":0.4267, "I-C8H18":0.3302, "C7H8":0.2431}.items()])
    X_fuel = "POSF10325:1.0"
    # create path to fuel model
    # fuel_model = os.path.join(DATA_DIR, "verification", "jet-fuel.yaml")
    # fuel_model = os.path.join(DATA_DIR, "farnesane", "farnesane.yaml")
    fuel_model = os.path.join(DATA_DIR, "A2NOx.yaml")
    # fuel_model = os.path.join(DATA_DIR, "combustor.yaml")
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
    return (combustor.thermo.state, thermo_states, x)


@click.command()
@click.option('--regenerate', is_flag=True, help='Regenerate verification data')
def combustor_verification(regenerate):
    convert_mission_out("CFM56_5B_EDB.out")
    # curve_fit_thrust_data(test=True)
    ver_dir = os.path.join(DATA_DIR, "verification")
    cbs = os.path.join(ver_dir, f"combustor-verification.hdf5")
    yfile = os.path.join(ver_dir, "combustor-verification.yaml")
    if not os.path.isdir(ver_dir):
        os.mkdir(ver_dir)
    fuel, equiv_ratio, X_fuel, X_air = combustor_design_params()
    thrust_levels = numpy.linspace(0.1, 1.0, 5)
    thrust_levels = [0.07, 0.3, 0.65, 0.85, 1.0]
    # thrust_levels = [1.0]
    if not os.path.isfile(cbs) or regenerate:
        # run in parallel
        with mp.Pool(min(os.cpu_count(), len(thrust_levels))) as pool:
            res = pool.map(parallel_run_combustor, thrust_levels)
        # solution array
        combustor_states = ct.SolutionArray(fuel, extra=["TL", "phi", "mdot_exhaust", "mdot_fuel"])
        verification_data = {}
        for r in res:
            if r is not None:
                state, thermo_states, tl = r
                verification_data[f"thrust-level-{tl}"] = thermo_states
                combustor_states.append(state, TL=tl, phi=equiv_ratio, mdot_fuel=thermo_states["mdot_fuel"], mdot_exhaust=thermo_states["mdot_out"])
        # dump verification data
        os.path.join("")
        with open(yfile, "w") as f:
            yaml.dump(verification_data, f)
        # write out state and reformats
        combustor_states.save(cbs, overwrite=True, name="thermo")
    # create plots of profiles from yaml
    with open(yfile, "r") as f:
        yaml_data = yaml.load(f, Loader=yaml.SafeLoader)
    # Mean equivalence ratio from cycle output
    fig, ax = plt.subplots()
    equiv_means = [get_prop_by_thrust_level(t, "FAR")/0.068 for t in thrust_levels]
    tls, phi = zip(*[(v["thrust_level"], v["phi_mean"]) for k, v in yaml_data.items()])
    ax.plot(tls, phi, color=COLORS[0])
    ax.plot(thrust_levels, equiv_means, linestyle="", marker="s", color=COLORS[0])
    # ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_ylabel("Mean Equivalence Ratio")
    ax.set_xlabel("Thrust Percentage")
    plt.grid(visible=True)
    fig.savefig(os.path.join(ver_dir, f"mean-equiv-ratio.pdf"), bbox_inches='tight')
    plt.close()
    # temperature profiles
    fig, ax = plt.subplots()
    for i, k in enumerate(yaml_data.keys()):
        ldata = yaml_data[k]
        ax.plot(ldata["phi_distribution"], ldata["temperature_distribution"], color=COLORS[i], label=f"{float(k.split('-')[-1]):.2f}", marker="o")
    ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_yticks(numpy.arange(500, 2500, 500))
    ax.set_ylabel("Temperature [K]")
    ax.set_xlabel("Equivalence Ratio")
    plt.grid(visible=True)
    fig.savefig(os.path.join(ver_dir, f"T-profiles.pdf"), bbox_inches='tight')
    plt.close()
    # mass flow fraction profiles
    with open (os.path.join(ver_dir, "mfeq-verification.csv")) as f:
        rdr = csv.reader(f)
        next(rdr)
        next(rdr)
        mfknown = [list(map(float, filter(lambda x: x, k))) for k in zip(*[r for r in rdr ])]
    # mass fraction verification
    fig, ax = plt.subplots()
    # plot all of the known lists
    for j, i in enumerate(range(0, len(mfknown), 2)):
        x = mfknown[i]
        y = mfknown[i+1]
        if i == 0:
            ax.plot(x, y, linestyle=":", color="k", marker="*", label="Source data")
        else:
            ax.plot(x, y, linestyle=":", color="k", marker="*")
    for i, k in enumerate(yaml_data.keys()):
        ldata = yaml_data[k]
        ax.plot(ldata["phi_distribution"], ldata["mass_flow_fraction_distribution"], color=COLORS[i], label=f"{float(k.split('-')[-1]):.2f}", marker="o")
    ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_ylabel("Fraction of Total Mass Flow")
    ax.set_xlabel("Equivalence Ratio")
    plt.grid(visible=True)
    ax.set_xlim([0, 3.5])
    ax.set_xticks(numpy.arange(0, 3.5, 0.5))
    fig.savefig(os.path.join(ver_dir, f"mass-flow-fraction-profiles.pdf"), bbox_inches='tight')
    # plt.show()
    plt.close()
    # create plot with generated data
    with h5py.File(cbs, "r") as hf:
        T_data = numpy.array(hf["thermo"]["data"]["T"])
        Y_data = hf["thermo"]["data"]["Y"]
        fuel_data = hf["thermo"]["data"]["mdot_fuel"]
        exhaust_data = hf["thermo"]["data"]["mdot_exhaust"]
        thrust_levels = numpy.array(hf["thermo"]["data"]["TL"])
        # compute NOx
        mass_nox = (Y_data[:, fuel.species_index("NO")] * exhaust_data) * 1000 # grams Y_data[:, fuel.species_index("NO2")] * exhaust_data
        nox_ratio = mass_nox / fuel_data
        mass_co = Y_data[:, fuel.species_index("CO")] * exhaust_data * 1000 # grams
        co_ratio = mass_co / fuel_data

    fig, axes = plt.subplots(1, 3)
    with open(os.path.join(ver_dir, "NOX-data.csv")) as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)
        nox_x, nox_y = zip(*list(map(lambda x: (float(x[0])/100, float(x[1])), reader)))
    # plot
    axes[0].plot(thrust_levels, nox_ratio, color=COLORS[0], label="Model prediction")
    axes[0].plot(nox_x, nox_y, color=COLORS[1], linestyle="", marker="o", label="EDB data")
    axes[0].set_xlabel("Thrust Fraction")
    axes[0].set_ylabel("EI $\\text{NO}_x$ [g/$\\text{kg}_{\\text{fuel}}$]")

    with open(os.path.join(ver_dir, "CO-data.csv")) as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)
        co_x, co_y = zip(*list(map(lambda x: (float(x[0])/100, float(x[1])), reader)))
    # plot
    axes[1].plot(thrust_levels, co_ratio, color=COLORS[0], label="Model prediction")
    axes[1].plot(co_x, co_y, color=COLORS[1], linestyle="", marker="o", label="EDB data")
    axes[1].set_xlabel("Thrust Fraction")
    axes[1].set_ylabel("EI CO [g/$\\text{kg}_{\\text{fuel}}$]")
    # plot
    axes[2].plot(thrust_levels, T_data, color=COLORS[0], label="Model prediction")
    axes[2].set_xlabel("Thrust Fraction")
    axes[2].set_ylabel("Temperature [K]")
    plt.subplots_adjust(wspace=0.25)
    fig.savefig(os.path.join(ver_dir, f"emissions-indices.pdf"), bbox_inches='tight')
    # plt.show()
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
    with open (os.path.join(ver_dir, "sz-temps-7.csv")) as f:
        rdr = csv.reader(f)
        next(rdr)
        next(rdr)
        zfvs, Tfvs, zsvs, Tsvs = [list(map(float, filter(lambda x: x, k))) for k in zip(*[r for r in rdr ])]
    # 100% thrust temperature profiles
    fig, ax = plt.subplots()
    ax.plot(zfv, Tfv, color=COLORS[0], marker="o", linestyle="", markerfacecolor="white")
    ax.plot(zsv, Tsv, color=COLORS[1], marker="o", linestyle="", markerfacecolor="white")
    ax.plot(z_fast_data, T_fast_data, color=COLORS[0])
    ax.plot(z_slow_data, T_slow_data, color=COLORS[1])
    plt.grid(visible=True)
    ax.set_xticks(numpy.arange(0, 120, 20))
    ax.set_xlim([-10, 110])
    ax.set_xlabel("Length [%]")
    ax.set_ylabel("Temperature [K]")
    custom_lines = [Line2D([0], [0], color=COLORS[0], lw=2),
                Line2D([0], [0], color=COLORS[1], lw=2),
                Line2D([0], [0], color="k", lw=0, marker="o"),]
    ax.legend(custom_lines, ["fast", "slow", "Source Data"], loc="upper right")
    fig.savefig(os.path.join(ver_dir, f"pfr-temperature-profile-100.pdf"), bbox_inches='tight')
    plt.close()
    # 7% thrust temperature profiles
    fig, ax = plt.subplots()
    ax.plot(z_fast_data_s, T_fast_data_s, color=COLORS[0])
    ax.plot(z_slow_data_s, T_slow_data_s, color=COLORS[1])
    ax.plot(zfvs, Tfvs, color=COLORS[0], marker="o", linestyle="")
    ax.plot(zsvs, Tsvs, color=COLORS[1], marker="o", linestyle="")
    plt.grid(visible=True)
    ax.set_xticks(numpy.arange(0, 100, 20))
    ax.set_xlim([0, 100.0])
    ax.set_xlabel("Length [%]")
    ax.set_ylabel("Temperature [K]")
    custom_lines = [Line2D([0], [0], color=COLORS[0], lw=2),
                Line2D([0], [0], color=COLORS[1], lw=2),
                Line2D([0], [0], color="k", lw=0, marker="o"),]
    ax.legend(custom_lines, ["fast", "slow", "Source Data"], loc="upper right")
    fig.savefig(os.path.join(ver_dir, f"pfr-temperature-profile-7.pdf"), bbox_inches='tight')
    plt.close()

def convert_mission_out(mout_name="CFM56_5B_737Mission.out"):
    ver_dir = os.path.join(DATA_DIR, "verification")
    fname = os.path.join(ver_dir, mout_name)
    with open(fname, "r") as f:
        lines = f.read().split("\n")
        headers = ",".join(re.split("\s+", lines[2])[1:-1])
        nls = [headers]
        for l in lines[3:]:
            r = re.split("\s+", l)[2:-1]
            nls.append(",".join(r))
        content = "\n".join(nls)
    with open(os.path.join(ver_dir, f"{mout_name.split('.')[0]}.csv"), "w") as f:
        f.write(content)

def check_thrust():
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", 'CFM56-7B-737.csv'), delimiter=",", engine="python")
    ma = cycle_data['W3[kg/s]'].values
    mf = cycle_data['Wf[kg/s]'].values
    thrust = cycle_data['NetThrust[kN]'].values
    for i, j, k in zip(ma, mf, thrust):
        print(i, j, k)

def initial_conditions_solver(Y):
    density = 1.2
    Y[0] = 30 * 48 / density / 1e6 # O3
    Y[1] = 100 * 28.01 / density / 1e6 # CO
    Y[2] = 1 - numpy.sum(Y[:2])
    return Y

def mcm_verification():
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    volume = 1e4
    gas.TPX = 298, ct.one_atm, "O2:0.205, N2:0.785, H2O:0.5"
    atm_air = ct.Reservoir(gas)
    Y = fsolve(initial_conditions_solver, [0.1, 0.1, 0.8])
    gas.Y = f"O3: {Y[0]}, CO:{Y[1]}, O2:{0.205*Y[2]}, N2:{0.785*Y[2]}, H2O:{0.01*Y[2]}"
    # setup initial conditions as 30 ppbv O3 and 100 ppbv CO with 1% water vapor
    r = PlumeReactor(gas, start_day=182, no_change_species=["H2O", "O2", "N2"])
    r.entrainment = False # Turn off entrainment
    r.volume = volume
    r.energy_enabled = False
    # mfc = ct.MassFlowController(atm_air, r, mdot=100)
    # setup reactor network
    net = ct.ReactorNet([r])
    net.max_steps = 1e8
    # The atmospheric integration struggles with the preconditioner.
    # Perhaps best to exclude preconditioning from the application for now
    # net.preconditioner = ct.AdaptivePreconditioner()
    net.derivative_settings = {"skip-flow-devices":True}
    # set tolerances
    net.rtol = 1e-12
    net.atol = 1e-20
    net.initialize()
    net.step()
    # setup time steps
    ndays = 3
    seconds = ndays * 24 * 3600 # convert days to seconds
    # setup solution array
    arr = ct.SolutionArray(gas, extra=["t", "mass", "volume"])
    arr.append(r.thermo.state, t=0.0, mass=r.mass, volume=r.volume)
    # loop through all time steps
    times = numpy.linspace(0, seconds, 5000)
    failed = False
    for i, t in enumerate(times[1:]): #for t in steps:
        failures = 0
        t_target = t
        net.max_time_step = 10
        while net.time < t:
            try:
                print(f"Advancing to {t:0.2f}..")
                former_time = net.time
                net.advance(t_target)
                # reset targeted time and time step
                t_target = t
                failures = 0
                net.max_time_step = 10
                # append the state so last state is always available
                arr.append(r.thermo.state, t=net.time, mass=r.mass, volume=r.volume)
            except Exception as e:
                print(f"Failure {failures}: reducing max time step")
                # reset reactor and network
                r.thermo.TDX = arr[-1].TDX
                r.syncState()
                r.volume = arr[-1].volume
                net.initial_time = former_time
                net.max_time_step = net.max_time_step / 10
                net.initialize()
                t_target = former_time + 1 # 10 seconds
                failures += 1
                if failures >= 20:
                    failed = True
                    break
        if failed:
            break
    # save hdf5 file
    hdf5_file = os.path.join(ver_dir, "mcm.hdf5")
    arr.save(hdf5_file, overwrite=True, name="thermo")
    reformat_hdf5(hdf5_file)
    # verification data
    mcm_data = pd.read_csv(os.path.join(ver_dir, 'mcm-verification.csv'), engine="python")
    mcm_time = mcm_data['TIME'].values
    # plot from arr object
    fig, axes = plt.subplots(2, 3)
    plt.subplots_adjust(wspace=1.25, hspace=0.5)
    # closure to repeat plots
    def plot_spec(ax, sp, ylabel="", scalar=1):
        mct = mcm_time / 24 / 3600
        mcs = 1
        mcmx = numpy.array([mct[i] for i in range(0, len(mcm_data[sp].values), mcs)])
        mcmy = numpy.array([mcm_data[sp].values[i] for i in range(0, len(mcm_data[sp].values), mcs)]);
        mcmy = mcmy / numpy.amax(mcmy)
        sp_id = gas.species_index(sp)
        model_data = arr.Y[:, sp_id] * arr.density / gas.molecular_weights[sp_id]
        model_data = model_data / numpy.amax(model_data)
        # plots
        ax.plot(mcmx, mcmy, color=COLORS[0], linestyle="", marker="o", markerfacecolor="white")
        ax.plot(arr.t / 24 / 3600, model_data , color=COLORS[1])
        ax.set_xlim(0, 3)
        ax.set_ylim(0, 1.2)
        ax.grid(True)
        ax.set_xticks(numpy.linspace(0, 1, 4))
        ax.set_xticks(numpy.linspace(0, 3, 4))
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(ylabel)
    plot_spec(axes[0][0], "O3", "$O_3$", arr.density / gas.molecular_weights[gas.species_index("O3")])
    plot_spec(axes[0][1], "CO", "$CO$", arr.density / gas.molecular_weights[gas.species_index("CO")])
    plot_spec(axes[0][2], "O1D", "$O_{1D}$", arr.density / gas.molecular_weights[gas.species_index("O1D")])
    plot_spec(axes[1][0], "OH", "$OH$", arr.density / gas.molecular_weights[gas.species_index("OH")])
    plot_spec(axes[1][1], "HO2", "$HO_2$", arr.density / gas.molecular_weights[gas.species_index("HO2")])
    plot_spec(axes[1][2], "H2O2", "$H_2O_2$", arr.density / gas.molecular_weights[gas.species_index("H2O2")])
    plt.savefig(os.path.join(ver_dir, "mcm-verification.pdf"))
    plt.close()

    # mcm temperature plot
    fig, ax = plt.subplots()
    ax.plot(mcm_time / 24 / 3600, mcm_data["TEMP"].values, color=COLORS[0], linestyle="", marker="o", markerfacecolor="white")
    ax.plot(arr.t / 24 / 3600, arr.T , color=COLORS[1])
    plt.savefig(os.path.join(ver_dir, "mcm-temp.pdf"))
    plt.close()

def test_solar_zenith_angle():
    # setup reactor
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    # setup initial conditions as 30 ppbv O3 and 100 ppbv CO with 1% water vapor
    r = PlumeReactor(gas, start_day=180)
    # one week
    hours = 168
    times = list(range(hours))
    angles = [numpy.rad2deg(r.calculate_solar_zenith_angle(t*3600)) for t in times] # a week of angles
    fig, ax = plt.subplots()
    stimes = [t*3600 for t in times]
    ax.plot(stimes, angles)
    ax.set_xticks([i * 3600 for i in range(0, 169, 24)])
    ax.set_xticklabels([str(i) for i in range(168//24 + 1)])
    plt.show()


if __name__ == "__main__":
    # get_mass_flow_thrust_data(0.5)
    check_thrust()
