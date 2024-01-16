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
from sklearn.metrics import mean_squared_error
from cac.constants import DATA_DIR, COLORS
from cac.combustor import multizone_combustor, get_prop_by_thrust_level, get_interpolated_property
from cac.combustor import curve_fit_thrust_data, reformat_hdf5
from cac.reactors import PlumeSolution, PlumeReactor
from cac.rates import ZenithAngleRate, ZenithAngleData

def rms_deviation_values(name, xF, yF, xV, yV):
    rms_file = os.path.join(DATA_DIR, "verification", "rmse.yaml")
    if os.path.isfile(rms_file):
        with open(rms_file, "r") as f:
            yaml_data = yaml.load(f, Loader=yaml.SafeLoader)
    else:
        yaml_data = {}
    yfunc = interp1d(xV, yV)
    keep = numpy.logical_and(xF>numpy.amin(xV), xF<numpy.amax(xV))
    points = list(filter(lambda x: keep[x], range(len(keep))))
    xF = xF[points[0]:points[-1]]
    yF = yF[points[0]:points[-1]]
    yvals = yfunc(xF)
    yaml_data[name] = float(mean_squared_error(yF, yvals, squared=False))
    with open(rms_file, "w") as f:
        yaml.dump(yaml_data, f)
    return yaml_data[name]

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
    plt.savefig(fname)
    # find rms error
    rms_deviation_values("box-moles", times, normalized_Xh2o, xX, yX)
    rms_deviation_values("box-temp", times, normalized_temp, xT, yT)

def ANOx_model():
    return "POSF10325:1.0", os.path.join(DATA_DIR, "A2NOx.yaml")

def creck_model():
    return "NC12H26:0.404, IC8H18:0.295, TMBENZ:0.073, NPBENZ:0.228", os.path.join(DATA_DIR, "verification", "jet-fuel.yaml")

def polimi_model():
    return "NC12H26:0.404, IC8H18:0.295, TMBENZ:0.073, NPBENZ:0.228", os.path.join(DATA_DIR, "verification", "polimi.yaml")

def hydogen_model():
    return "H2:1.0", "h2o2.yaml"

def combustor_design_params():
    # design parameters
    equiv_ratio = 1.77
    EPR = 32.6 # engine pressure ratio
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    X_air = "O2:0.21, N2:0.79"
    # create path to fuel model
    X_fuel, fuel_model = polimi_model()
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
    out_dir = os.path.join(DATA_DIR, "verification")
    try:
        print(f"Running thrust-level: {x:.3f}")
        combustor, thermo_states = multizone_combustor(fuel, x, equiv_ratio, X_fuel, X_air, n_pz=21, name_id=f"-verification-{x:0.2f}", outdir=out_dir)
    except Exception as e:
        print(f"Failed for Thrust-level: {x:.3f}")
        print(e)
        return None
    return (combustor.thermo.state, thermo_states, x)

def generate_mixing_zone_plot(thrust_level, source_file=None, variable="T", label="Temperature [K]"):
    ver_dir = os.path.join(DATA_DIR, "verification")
    # create figure and plot
    fig, ax = plt.subplots()
    # add source data plot
    if source_file is not None:
        with open (os.path.join(ver_dir, source_file)) as f:
            rdr = csv.reader(f)
            next(rdr)
            next(rdr)
            zfv, Tfv, zsv, Tsv = [list(map(float, filter(lambda x: x, k))) for k in zip(*[r for r in rdr ])]
            ax.plot(zfv, Tfv, color=COLORS[0], marker="o", linestyle="", markerfacecolor="white")
            ax.plot(zsv, Tsv, color=COLORS[1], marker="o", linestyle="", markerfacecolor="white")
    # model prediction data
    # fast mixing data
    fmh5 = os.path.join(ver_dir, f"fast-mixing-states-verification-{thrust_level:.2f}.hdf5")
    with h5py.File(fmh5, "r") as fhf:
        z_fast_data = numpy.array(fhf["z"])
        z_fast_data /= numpy.amax(z_fast_data)
        z_fast_data *= 100
        T_fast_data = numpy.array(fhf[variable])
    # slow mixing data
    smh5 = os.path.join(ver_dir, f"slow-mixing-states-verification-{thrust_level:.2f}.hdf5")
    with h5py.File(smh5, "r") as shf:
        z_slow_data = numpy.array(shf["z"])
        z_slow_data /= numpy.amax(z_slow_data)
        z_slow_data *= 100
        T_slow_data = numpy.array(shf[variable])
    # add model data
    ax.plot(z_fast_data, T_fast_data, color=COLORS[0])
    ax.plot(z_slow_data, T_slow_data, color=COLORS[1])
    plt.grid(visible=True)
    ax.set_xticks(numpy.arange(0, 120, 20))
    ax.set_xlim([-10, 110])
    ax.set_xlabel("Length [%]")
    ax.set_ylabel(label)
    custom_lines = [Line2D([0], [0], color=COLORS[0], lw=2),
                Line2D([0], [0], color=COLORS[1], lw=2)]
    labels = ["Fast Mixing", "Slow Mixing"]
    if source_file is not None:
        custom_lines.append(Line2D([0], [0], color="k", lw=0, marker="o"))
        labels.append("Source Data")
    ax.legend(custom_lines, labels, loc="upper right")
    fig.savefig(os.path.join(ver_dir, f"pfr-{variable}-profile-{thrust_level*100:.0f}.pdf"), bbox_inches='tight')
    plt.close()

def generate_ei_chart(EI, fuel, cbs):
    # create plot with generated data
    ver_dir = os.path.join(DATA_DIR, "verification")
    with h5py.File(cbs, "r") as hf:
        Y_data = hf["thermo"]["data"]["Y"]
        fuel_data = hf["thermo"]["data"]["mdot_fuel"]
        exhaust_data = hf["thermo"]["data"]["mdot_exhaust"]
        thrust_levels = numpy.array(hf["thermo"]["data"]["TL"])
        # compute EI
        if EI == "NOx":
            mass_ei = (Y_data[:, fuel.species_index("NO")] + Y_data[:, fuel.species_index("NO2")]) * exhaust_data * 1000
            ei_ratio = mass_ei / fuel_data
        else:
            mass_ei = (Y_data[:, fuel.species_index(EI)]) * exhaust_data * 1000
            ei_ratio = mass_ei / fuel_data
    # plot
    fig, ax = plt.subplots()
    ax.plot(thrust_levels, ei_ratio, color=COLORS[0], label="Model prediction")
    # plot experimental data
    icao_data = pd.read_csv(os.path.join(DATA_DIR, "verification", 'icao-verify-7B.csv'), delimiter=",", engine="python")
    max_thrust = get_prop_by_thrust_level(1.0, 'NetThrust[kN]')
    thrusts = numpy.array([get_interpolated_property(v, 'Wf[kg/s]', 'NetThrust[kN]') for v in icao_data["mdotf"].values]) / max_thrust
    # plot
    ax.plot(thrusts, icao_data[EI].values, color=COLORS[1], linestyle="", marker="o", label="EDB data", markerfacecolor="white")
    ax.set_xlabel("Thrust Fraction")
    ax.set_ylabel(f"EI {EI} "+" [g/$\\text{kg}_{\\text{fuel}}$]")
    fig.savefig(os.path.join(ver_dir, f"EI-{EI}.pdf"), bbox_inches='tight')
    plt.close()

def generate_pz_phi_plot(variable, label, add_source=False):
    ver_dir = os.path.join(DATA_DIR, "verification")
    yfile = os.path.join(ver_dir, "combustor-verification.yaml")
    with open(yfile, "r") as f:
        yaml_data = yaml.load(f, Loader=yaml.SafeLoader)
    if add_source:
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
                ax.plot(x, y, linestyle=":", color="k", marker="o", label="Source data", markerfacecolor="white")
            else:
                ax.plot(x, y, linestyle=":", color="k", marker="o", markerfacecolor="white")
    else:
        fig, ax = plt.subplots()
    for i, k in enumerate(yaml_data.keys()):
        ldata = yaml_data[k]
        ax.plot(ldata["phi_distribution"], ldata[variable], color=COLORS[i], label=f"{float(k.split('-')[-1]):.2f}", marker="o")
    ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_ylabel(label)
    ax.set_xlabel("Equivalence Ratio")
    plt.grid(visible=True)
    save_label = variable.replace("_", "-")
    fig.savefig(os.path.join(ver_dir, f"{save_label}-profiles.pdf"), bbox_inches='tight')
    plt.close()

def generate_mean_equiv_ratio_plot():
    ver_dir = os.path.join(DATA_DIR, "verification")
    yfile = os.path.join(ver_dir, "combustor-verification.yaml")
    with open(yfile, "r") as f:
        yaml_data = yaml.load(f, Loader=yaml.SafeLoader)
    # Mean equivalence ratio from cycle output
    fig, ax = plt.subplots()
    tls, phi = zip(*[(v["thrust_level"], v["phi_mean"]) for k, v in yaml_data.items()])
    equiv_means = [get_prop_by_thrust_level(t, "FAR")/0.068 for t in tls]
    ax.plot(tls, phi, color=COLORS[0])
    ax.plot(tls, equiv_means, linestyle="", marker="s", color=COLORS[0])
    # ax.legend(ncols=1, bbox_to_anchor=(1.2, 1.03))
    ax.set_ylabel("Mean Equivalence Ratio")
    ax.set_xlabel("Thrust Percentage")
    plt.grid(visible=True)
    fig.savefig(os.path.join(ver_dir, f"mean-equiv-ratio.pdf"), bbox_inches='tight')
    plt.close()

@click.command()
@click.option('--regenerate', is_flag=True, help='Regenerate verification data')
def combustor_verification(regenerate):
    ver_dir = os.path.join(DATA_DIR, "verification")
    # set model function
    cbs = os.path.join(ver_dir, f"combustor-verification.hdf5")
    yfile = os.path.join(ver_dir, "combustor-verification.yaml")
    if not os.path.isdir(ver_dir):
        os.mkdir(ver_dir)
    fuel, equiv_ratio, X_fuel, X_air = combustor_design_params()
    thrust_levels = numpy.linspace(0.07, 1.0, 20)
    thrust_levels = [0.07, 0.3, 0.65, 0.85, 1.0]
    # thrust_levels = [0.07, 1.0]
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
    # functions to generate plots
    generate_mean_equiv_ratio_plot()
    # primary zone EI
    generate_pz_phi_plot("mass_flow_fraction_distribution", "Mass Flow Fraction", add_source=True)
    generate_pz_phi_plot("temperature_distribution", "Temperature [K]")
    generate_pz_phi_plot("EI_CO_pz", "EI CO $g / kg_{fuel}$")
    generate_pz_phi_plot("EI_NOx_pz", "EI NOx $g / kg_{fuel}$")
    # Total emissions indices
    generate_ei_chart("NOx", fuel, cbs)
    generate_ei_chart("CO", fuel, cbs)
    # mixing plots
    generate_mixing_zone_plot(1.0, "sz-T-100.csv")
    generate_mixing_zone_plot(0.07, "sz-T-7.csv")
    generate_mixing_zone_plot(1.0, variable="phi", label="Equivalence Ratio")
    generate_mixing_zone_plot(0.07, variable="phi", label="Equivalence Ratio")

def convert_mission_out(mout_name="CFM56_7B_737Mission.out"):
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
    Y[0] = 30 * 48 / 1e9 * density# O3
    Y[1] = 100 * 28.01 / 1e9 * density # CO
    remainder = 1 - numpy.sum(Y[:2])
    Y[2] = 0.205 * remainder
    Y[3] = 0.785 * remainder
    Y[4] = 0.01 * remainder
    Y[5] = 1 - numpy.sum(Y)
    return Y

def mcm_verification():
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    volume = 1e4
    gas.TPX = 298, ct.one_atm, "O2:0.205, N2:0.785, H2O:0.5"
    atm_air = ct.Reservoir(gas)
    Y = fsolve(initial_conditions_solver, [0.1, 0.1, 0.2, 0.7, 0.01, 1])
    gas.Y = f"O3: {Y[0]}, CO:{Y[1]}, O2:{Y[2]}, N2:{Y[3]}, H2O:{Y[4]}"
    # setup initial conditions as 30 ppbv O3 and 100 ppbv CO with 1% water vapor
    r = PlumeReactor(gas, start_day=182, altitude=1e3, no_change_species=["H2O", "O2", "N2"])
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
    rates = numpy.hstack([[0], r.kinetics.forward_rate_constants])
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
                trts = numpy.hstack([[net.time], r.kinetics.forward_rate_constants])
                rates = numpy.vstack([rates, trts])
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
    csv_file = os.path.join(ver_dir, "mcm.csv")
    arr.save(hdf5_file, overwrite=True, name="thermo")
    arr.save(csv_file, overwrite=True)
    reformat_hdf5(hdf5_file)
    # save rates file
    rates_file = os.path.join(ver_dir, "rates.csv")
    with open(rates_file, "w") as f:
        cwriter = csv.writer(f, delimiter=",")
        cwriter.writerow(["time"]+r.kinetics.reactions())
        for rt in rates:
            cwriter.writerow(rt)
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
        ax.plot(mcmx, mcmy, color=COLORS[0], linestyle="", marker="o", markerfacecolor="white", label="MCM")
        ax.plot(arr.t / 24 / 3600, model_data , color=COLORS[1], label="Cantera")
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
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncols=2)
    plt.savefig(os.path.join(ver_dir, "mcm-verification.pdf"))
    plt.close()

    # mcm temperature plot
    fig, ax = plt.subplots()
    ax.plot(mcm_time / 24 / 3600, mcm_data["TEMP"].values, color=COLORS[0], linestyle="", marker="o", markerfacecolor="white")
    ax.plot(arr.t / 24 / 3600, arr.T , color=COLORS[1])
    plt.savefig(os.path.join(ver_dir, "mcm-temp.pdf"))
    plt.close()

    # RMSE
    for sp in ["O3", "O1D", "OH", "H2O2", "HO2", "CO"]:
        mcmy = mcm_data[sp].values
        mcmy = mcmy / numpy.amax(mcmy)
        sp_id = gas.species_index(sp)
        model_data = arr.Y[:, sp_id] * arr.density / gas.molecular_weights[sp_id]
        model_data = model_data / numpy.amax(model_data)
        rms_deviation_values(f"mcm-{sp}", arr.t, model_data, mcm_time, mcmy)

def test_solar_zenith_angle():
    # setup reactor
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    sd = 182
    # setup initial conditions as 30 ppbv O3 and 100 ppbv CO with 1% water vapor
    r = PlumeReactor(gas, start_day=sd, altitude=1e3)
    # one week
    ndays = 365
    hours = ndays * 24
    times = list(range(hours))
    angles = [numpy.rad2deg(r.calculate_solar_zenith_angle(t*3600)) for t in times] # a week of angles
    fig, ax = plt.subplots()
    stimes = [t / 24 + sd for t in times]
    ax.plot(stimes, angles)
    ax.set_ylim([-10, 190])
    ax.set_yticks(numpy.arange(0, 190, 20))
    days = numpy.linspace(0+sd, ndays+sd, 10)
    days = [int(day) for day in days]
    ax.set_xticks(days)
    plt.savefig(os.path.join(ver_dir, "zenith-angles.pdf"))
    plt.close()
    # plot solar zenith maximums
    fig, ax = plt.subplots()
    times = numpy.arange(12, hours, 24)
    angles = [numpy.rad2deg(r.calculate_solar_zenith_angle(t*3600)) for t in times]
    stimes = [t / 24 + sd for t in times]
    ax.plot(stimes, angles)
    ax.set_ylim([-10, 190])
    ax.set_yticks(numpy.arange(0, 190, 20))
    days = numpy.linspace(0+sd, ndays+sd, 10)
    days = [int(day) for day in days]
    ax.set_xticks(days)
    plt.savefig(os.path.join(ver_dir, "zenith-angles-maximums.pdf"))
    plt.close()
    # plot solar zenith minimums
    fig, ax = plt.subplots()
    times = numpy.arange(0, hours, 24)
    angles = [numpy.rad2deg(r.calculate_solar_zenith_angle(t*3600)) for t in times]
    stimes = [t / 24 + sd for t in times]
    ax.plot(stimes, angles)
    ax.set_ylim([-10, 190])
    ax.set_yticks(numpy.arange(0, 190, 20))
    days = numpy.linspace(0+sd, ndays+sd, 10)
    days = [int(day) for day in days]
    ax.set_xticks(days)
    plt.savefig(os.path.join(ver_dir, "zenith-angles-minimums.pdf"))
    plt.close()
    # output of J1OD and JNO2 at local solar noon
    photo_rates = {}
    gas.zenith_angle = r.calculate_solar_zenith_angle(12*3600)
    # custom zenith data
    zdata = ZenithAngleData()
    zdata.update(gas)
    # setup yaml
    photo_rates["day"] = sd
    photo_rates["altitude"] = float(r.altitude)
    photo_rates["zenith-angle"] = float(numpy.rad2deg(r.thermo.zenith_angle))
    photo_rates["rates"] = {}
    # loop over photo reactions
    for i, cr in enumerate(r.kinetics.reactions()):
        if isinstance(cr.rate, ZenithAngleRate):
            photo_rates["rates"][str(cr)] = {}
            photo_rates["rates"][str(cr)]["J"] = float(cr.rate.eval(zdata))
            params = {}
            cr.rate.get_parameters(params)
            params = {k:float(p) for k,p in params.items()}
            photo_rates["rates"][str(cr)].update(params)
    with open(os.path.join(ver_dir, "J-rates.yaml"), "w") as f:
        yaml.dump(photo_rates, f)

def mcm_arrhenius_rate():
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    volume = 1e4
    T1 = 298
    P1 = ct.one_atm
    gas.TPX = T1, P1, "O2:0.205, N2:0.785, H2O:0.5"
    atm_air = ct.Reservoir(gas)
    Y = fsolve(initial_conditions_solver, [0.1, 0.1, 0.2, 0.7, 0.01, 1])
    gas.Y = f"O3: {Y[0]}, CO:{Y[1]}, O2:{Y[2]}, N2:{Y[3]}, H2O:{Y[4]}"
    r = PlumeReactor(gas, start_day=182, altitude=1e3, no_change_species=["H2O", "O2", "N2"])
    r.entrainment = False # Turn off entrainment
    r.volume = volume
    r.energy_enabled = False
    # setup reactor network
    net = ct.ReactorNet([r])
    net.advance(12*2600)
    # Reaction 15: O + O3 => 2 O2
    for i, rt in enumerate(r.kinetics.reactions()):
        if rt.__str__() == "O + O3 => 2 O2":
            reaction = rt
            rid = i
    Mvalue = r.pressure() / (r.T * ct.boltzmann) / (100**3)
    cantera_rate = r.kinetics.forward_rate_constants[rid] * 100 ** 3 / ct.avogadro
    # calculated rate - 8.0D-12*EXP(-2060/TEMP) from fac file
    calc_rate = 8e-12 * numpy.exp(-2060 / T1)
    # cantera_rate = reaction.rate
    assert numpy.isclose(calc_rate, cantera_rate)

def mcm_kunknown_rate():
    ver_dir = os.path.join(DATA_DIR, "verification")
    gas = PlumeSolution(os.path.join(ver_dir, "verification.yaml"))
    volume = 1e4
    T1 = 298
    P1 = ct.one_atm
    gas.TPX = T1, P1, "O2:0.205, N2:0.785, H2O:0.5"
    atm_air = ct.Reservoir(gas)
    Y = fsolve(initial_conditions_solver, [0.1, 0.1, 0.2, 0.7, 0.01, 1])
    gas.Y = f"O3: {Y[0]}, CO:{Y[1]}, O2:{Y[2]}, N2:{Y[3]}, H2O:{Y[4]}"
    r = PlumeReactor(gas, start_day=182, altitude=1e3, no_change_species=["H2O", "O2", "N2"])
    r.entrainment = False # Turn off entrainment
    r.volume = volume
    r.energy_enabled = False
    # setup reactor network
    net = ct.ReactorNet([r])
    net.advance(12*2600)
    # both reactions of interest
    reactions = list(filter(lambda x: x[1].__str__() == "O + O2 => O3", enumerate(r.kinetics.reactions())))
    for er in reactions:
        print(er)
    conv = 100 ** 3 / ct.avogadro
    kfs = [r.kinetics.forward_rate_constants[i] * conv for i, rt in reactions]
    Mvalue = r.pressure() / (r.T * ct.boltzmann) / (100**3)
    N2 = 0.785 * Mvalue
    O2 = 0.205 * Mvalue
    # calculated rates
    # % 5.6D-34*N2*(TEMP/300)@(-2.6)*O2 : O = O3 ;
    # % 6.0D-34*O2*(TEMP/300)@(-2.6)*O2 : O = O3 ;
    calc_rates = [5.6e-34 * N2 * (T1 / 300) ** (-2.6) * O2, 6.0e-34 * O2 * (T1 / 300) ** (-2.6) * O2]
    calc_rates.sort()
    kfs.sort()
    assert numpy.allclose(kfs, calc_rates, rtol=1e-4)
    print(calc_rates)
    print(kfs)
    print(gas.concentrations[gas.species_index("O2")])

def function_tester():
    print("Replace me with whatever function you want to test and run verify_tester")
    mcm_kunknown_rate()
