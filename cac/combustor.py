import os
import sys
import csv
import h5py
import yaml
import click
import numpy
import warnings
import pandas as pd
import cantera as ct
import scipy.stats as stats
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from cac.merger import map_models_wrapped
from cac.constants import DATA_DIR
from cac.reactors import PlumeSolution, PlumeReactor, DilutionReactor

PRECONDITIONED = True
HIGH_TOLERANCE = True
ENGINE_MODEL = "7B"
def reformat_hdf5(hf_name):
    with h5py.File(hf_name, "r+") as hf:
        arr = hf["thermo"]['data']
        Y_data = hf["thermo"]["data"]["Y"]
        nsp = Y_data.shape[1]
        Y_arr = Y_data[()]
        # create species data sets
        for i, c in enumerate(arr.attrs["components"][2:2+nsp]):
            try:
                del hf[f"Y_{c}"] # delete it if it exists
            except:
                pass
            hf.create_dataset(f"Y_{c}", data=Y_arr[:, i])
        # write the rest of the data out
        other_keys = list(arr.keys())
        other_keys.remove("Y")
        for n in other_keys:
            data = hf["thermo"]["data"][n][()]
            try:
                del hf[n] # delete it if it exists
            except:
                pass
            hf.create_dataset(n, data=data)
        # delete thermo group
        del hf["thermo"]

def get_debug_atmosphere():
    spec = ct.Species.list_from_file(os.path.join(DATA_DIR, "atmosphere.yaml"))
    spec_gas = PlumeSolution(thermo='ideal-gas', kinetics="gas", species=spec)
    ys = """
    - A: 7.23e+08
      Ea: 0.0
      b: 0.0
      equation: SO3 + H2O => SA
      pyfile: atmosphere_complex_rates.py
      ro2file: atmosphere-ro2-sum.txt
      species-names: H2O
      type: complex-rate
      """
    R = ct.Reaction.list_from_yaml(ys, spec_gas)
    for r in R:
        spec_gas.add_reaction(r)
    atmosphere = PlumeReactor(spec_gas)
    return atmosphere


def pressure(r):
    return r.T * ct.gas_constant * r.mass / r.thermo.mean_molecular_weight / r.volume


def print_reactor_stats(net, r):
        print(f"Integrated to {net.time}..")
        print(r.mass)
        print(r.T)
        print(r.volume)
        print(r.density)
        print(pressure(r))
        print()


def new_network(r, precon=None):
    if precon is None:
        precon = PRECONDITIONED
    # setup reactor network
    net = ct.ReactorNet(r)
    net.derivative_settings = {"skip-falloff": False,
                            "skip-third-bodies": False,
                            "skip-coverage-dependence": True,
                            "skip-electrochemistry": True,
                            "skip-flow-devices": True,
                            "skip-walls": True}
    net.max_steps = 5e6
    if HIGH_TOLERANCE:
        net.rtol = 1e-16
        net.atol = 1e-20
    # setup preconditioner
    if precon:
        net.preconditioner = ct.AdaptivePreconditioner()
    net.initialize()
    return net


def print_state_TP(T, P, i):
    print()
    print(80 * "*")
    tp_str = f"STATE {i}: {T:.2f}, {P:.0f}"
    ast_len = (80 - len(tp_str)) // 2
    print(ast_len * "*" + tp_str + ast_len * "*")
    print(80 * "*")
    print()


def steady_state_advance(net, reactor, solarr):
        # get default tolerances:
        atol = net.atol
        residual_threshold = 1e-4
        max_steps = 10000
        max_state_values = net.get_state()  # denominator for feature scaling
        for step in range(max_steps):
            previous_state = net.get_state()
            solarr.append(reactor.thermo.state, U=0.0, z=0.0, t=net.time)
            # take 10 steps (just to increase speed)
            for n1 in range(10):
                net.step()
            state = net.get_state()
            max_state_values = numpy.maximum(max_state_values, state)
            # determine feature_scaled residual
            residual = numpy.linalg.norm((state - previous_state)
                / (max_state_values + atol)) / numpy.sqrt(net.n_vars)
            if residual < residual_threshold:
                break
        if step == max_steps - 1:
            raise ct.CanteraError('Maximum number of steps reached before'
                               ' convergence below maximum residual')


def adjust_volume_to_preserve_mass(r, tm):
    # adjust volume such that mass is
    step = 1e-6
    vp = False
    vm = False
    while not numpy.isclose(r.mass, tm, rtol=1e-12, atol=1e-16):
        if r.mass > tm:
            r.volume -= step
            vm = True
        elif r.mass < tm:
            r.volume += step
            vp = True
        # adjust step if both were visited
        if vm and vp:
            step /= 2
            vm = False
            vp = False


def get_prop_by_thrust_level(thrust_level, prop):
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", f'CFM56-{ENGINE_MODEL}-737.csv'), delimiter=",", engine="python")
    propty = cycle_data[prop].values
    thrust = cycle_data['NetThrust[kN]']
    thrust_percents = thrust / numpy.amax(thrust)
    zipped = list(zip(thrust_percents, propty))
    thrust_percents, propty = zip(*sorted(zipped, key=lambda x: x[0]))
    for i in range(len(thrust_percents)):
        if numpy.isclose(thrust_level, thrust_percents[i], rtol=1e-2):
            return propty[i]
    # interpolate for values
    propty_func = interp1d(thrust_percents, propty)
    return propty_func(thrust_level)


def get_mass_flow_thrust_data(thrust_level):
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", f'CFM56-{ENGINE_MODEL}-737.csv'), delimiter=",", engine="python")
    ma = cycle_data['W3[kg/s]'].values
    mf = cycle_data['Wf[kg/s]'].values
    thrust = cycle_data['NetThrust[kN]']
    thrust_percents = thrust / numpy.amax(thrust)
    zipped = list(zip(thrust_percents, thrust, mf, ma))
    thrust_percents, thrust, mf, ma = zip(*sorted(zipped, key=lambda x: x[0]))
    for i in range(len(thrust_percents)):
        if numpy.isclose(thrust_level, thrust_percents[i], rtol=1e-2):
            return mf[i], ma[i], mf[-1], ma[-1]
    # interpolate for values
    mdot_fuel_func = interp1d(thrust_percents, mf)
    mdot_air_func = interp1d(thrust_percents, ma)
    thrust_func = interp1d(thrust_percents, thrust)
    thrust = thrust_func(thrust_level)
    mdot_fuel = mdot_fuel_func(thrust_level)
    mdot_air = mdot_air_func(thrust_level)
    mdot_fuel_to = mf[-1]
    mdot_air_to = ma[-1]
    return mdot_fuel, mdot_air, mdot_fuel_to, mdot_air_to


def get_interpolated_property(v1, prop1, prop2):
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", f'CFM56-{ENGINE_MODEL}-737.csv'), delimiter=",", engine="python")
    p1 = cycle_data[prop1].values
    p2 = cycle_data[prop2].values
    zipped = list(zip(p1, p2))
    p1, p2 = zip(*sorted(zipped, key=lambda x: x[0]))
    for i in range(len(p1)):
        if numpy.isclose(v1, p1[i], rtol=1e-2):
            return p2[i]
    interp_func = interp1d(p1, p2)
    return interp_func(v1)


def get_mean_equivalence_ratio_from_data(thrust_level):
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", 'equiv_mean.csv'), delimiter=",", engine="python")
    phi_values = cycle_data["phi_mean"].values
    thrust_values = cycle_data["thrust_level"].values
    zipped = list(zip(thrust_values, phi_values))
    thrust_values, phi_values = zip(*sorted(zipped, key=lambda x: x[0]))
    for i in range(len(thrust_values)):
        if numpy.isclose(thrust_level, thrust_values[i], rtol=1e-2):
            return phi_values[i]
    interp_func = interp1d(thrust_values, phi_values)
    return interp_func(thrust_level)


def multizone_combustor(fuel, thrust_level, equiv_pz, X_fuel, X_air, **kwargs):
    # default parameters
    thermo_states = kwargs.get("thermo_states", {"T":[], "P":[]})
    n_pz = kwargs.get("n_pz", 21)
    n_pz += 1 if n_pz % 2 == 0  else 0 # ensure always odd
    mixing_param = kwargs.get("mixing_param", 0.55) # mixing parameter 0.3
    volume = kwargs.get("volume", 0.05) # m^3 0.043
    mdot_fuel, mdot_air, mdot_fuel_takeoff, mdot_air_takeoff = get_mass_flow_thrust_data(thrust_level)
    # initial pressure and temperature
    P_i = get_prop_by_thrust_level(thrust_level, "Pt3[kPa]") * 1000
    T_i = get_prop_by_thrust_level(thrust_level, "Tt3[K]")
    print_state_TP(T_i, P_i, f"cst:{thrust_level:0.2f}")
    thermo_states["T"] = thermo_states.get("T", []) + [f"{T_i:0.2f}"]
    thermo_states["P"] = thermo_states.get("P", []) + [f"{P_i:0.2f}"]
    assert thrust_level <= 1.0
    equiv_sec_zone = kwargs.get("equiv_sz", 0.4 * equiv_pz)
    outdir = kwargs.get("outdir", "./")
    name_id = kwargs.get("name_id", "")
    lhv_data = kwargs.get("lhv_data", 43.5e6) # Lower heating value of fuel in data set J/kg
    lhv_model = kwargs.get("lhv_model", 43.1e6)
    fuel_scalar = kwargs.get("fuel_scalar", lhv_data / lhv_model)
    mdot_fuel *= fuel_scalar
    # fixed params from study
    A_sz = kwargs.pop("A_sz", 0.15) # meters squared
    L_sz = 0.075 # m
    f_sm = 0.5 # fraction of mass flow in slow mixing
    f_fm = 1 - f_sm
    l_sa_sm = 0.94 # 0.75
    l_sa_fm = 0.015 # 0.09
    l_dae = 1
    l_das = 0.95
    # different fueling options for equiv ratio
    if (kwargs.get("use_total_equiv", False)):
        print("Using total equivalence ratio option")
        phi_across = equiv_pz
        f_air_pz = {"7B":0.256, "5B":0.256}.get(ENGINE_MODEL, 0.256)
        f_air_sa = {"7B":0.635, "5B":0.635}.get(ENGINE_MODEL, 0.635)
        f_air_da = (1 - f_air_sa - f_air_pz)
        # keys for mixture
        akeys = [fuel.species_index(x.split(":")[0].strip()) for x in X_air.split(",")]
        fkeys = [fuel.species_index(x.split(":")[0].strip()) for x in X_fuel.split(",")]
        # stoichiometric fuel/air
        fuel.TP = T_i, P_i
        fuel.set_equivalence_ratio(1.0, X_fuel, X_air)
        mdot_st = numpy.sum(mdot_fuel / fuel.Y[fkeys])
        mdot_air_st= numpy.sum(mdot_st * fuel.Y[akeys])
        FARst = mdot_fuel / mdot_air_st
        # total air flow
        mdot_air =  mdot_air_st / phi_across
        FAR = mdot_fuel / mdot_air
        # find primary zone equivalence ratio
        equiv_pz = mdot_fuel / (mdot_air * f_air_pz) / FARst
        equiv_sec_zone = 0.4 * equiv_pz
    else:
        FARst = kwargs.get("FARst", 0.068) # This value substantially impacts EIs
        # adjustments to airflow can better fit thesis data but overall behavior seems
        # to mostly correct
        FAR = get_prop_by_thrust_level(thrust_level, "FAR")
        # fraction of total air
        f_air_pz = mdot_fuel_takeoff * fuel_scalar / (mdot_air_takeoff * equiv_pz * FARst)
        f_air_sa = mdot_fuel_takeoff / (equiv_sec_zone * FARst * mdot_air_takeoff)
        f_air_da = (1 - f_air_sa - f_air_pz)
    # check the phi across
    phi_across = FAR / FARst
    # check the dilution air
    if f_air_da < 0:
        print("Not enough dilution air, f_air_da set to 0")
        f_air_da = 0
    # update mass flow rates
    mdot = mdot_air * f_air_pz + mdot_fuel
    # now determine mean equiv ratio with takeoff air and mdot fuel
    equiv_mean = mdot_fuel / (mdot_air * f_air_pz * FARst)
    # calculated parameters
    beta_sa_sm = f_air_sa * f_sm * mdot_air / (l_sa_sm * L_sz)
    beta_sa_fm = f_air_sa * f_fm * mdot_air / (l_sa_fm * L_sz)
    beta_da = f_air_da * mdot_air / (l_dae - l_das) / L_sz
    # calculate phi values
    equiv_min = kwargs.get("equiv_min", 0.3)
    equiv_max = kwargs.get("equiv_max", 4.5)
    sigma = mixing_param * equiv_mean
    coeff = min([abs(equiv_mean - equiv_min) / sigma, abs(equiv_max - equiv_mean) / sigma, 2])
    if n_pz > 1:
        phis = numpy.linspace(equiv_mean - coeff * sigma, equiv_mean + coeff * sigma, n_pz)
    else:
        phis = numpy.array([equiv_mean])
    # calculate mass flows
    mass_flow_fractions = (1 / (sigma * numpy.sqrt(2*numpy.pi))) * numpy.exp(- (phis - equiv_mean) ** 2 / (2 * sigma ** 2 ))
    if len(phis) > 1:
        mass_flow_fractions *= (phis[1] - phis[0])
    else:
        mass_flow_fractions *= phis[0]
    # scale mass flow fractions to equal all of mdot
    mass_flow_fractions *= 1 / numpy.sum(mass_flow_fractions)
    mfs = mass_flow_fractions * mdot
    assert numpy.isclose(numpy.sum(mfs), mdot)
    # add to thermo states
    split_volumes = volume * mass_flow_fractions # split over zones
    # add all calculated parameters to output
    thermo_states["FAR"] = float(FAR)
    thermo_states["FARst"] =  float(FARst)
    thermo_states["mdot_fuel"] = float(mdot_fuel)
    thermo_states["mdot_air"] = float(mdot_air)
    thermo_states["mdot_in"] = float(mdot)
    thermo_states["f_air_pz"] = float(f_air_pz)
    thermo_states["f_air_sa"] = float(f_air_sa)
    thermo_states["use_total_equiv"] = kwargs.get("use_total_equiv", False)
    thermo_states["phi_sz"] = float(equiv_sec_zone)
    thermo_states["phi_pz"] = float(equiv_pz)
    thermo_states["phi_mean"] = float(equiv_mean)
    thermo_states["phi_across"] = float(phi_across)
    thermo_states["thrust_level"] = float(thrust_level)
    thermo_states["phi_distribution"] = [float(p) for p in phis]
    thermo_states["mass_flow_fraction_distribution"] = [float(p) for p in mass_flow_fractions]
    thermo_states["mass_flow_distribution"] = [float(p) for p in mfs]
    thermo_states["n-deviations"] = float(coeff)
    thermo_states["n-reactors"] = int(n_pz)
    # for testing purposes
    if kwargs.get("norun", False):
        # output the initial composition
        fuel.TP = T_i, P_i
        fuel.set_equivalence_ratio(equiv_mean, X_fuel, X_air, basis="mole")
        # attempt to make directory
        try:
            os.mkdir("initials")
        except:
            pass
        initial_states = ct.SolutionArray(fuel, extra=["time", "moles", "mass", "volume"])
        nmoles = P_i * volume / ct.gas_constant / T_i
        mass = nmoles * fuel.mean_molecular_weight
        initial_states.append(fuel.state, time=0, volume=volume, moles=nmoles, mass=mass)
        # write out state and reformat
        is5 = os.path.join("initials", f"initial-states-{name_id}.hdf5")
        initial_states.save(is5, overwrite=True, name="thermo")
        reformat_hdf5(is5)
        states_file = os.path.join("initials", f"thermo-states-{name_id}.yaml")
        with open(states_file, "w") as f:
            yaml.dump(thermo_states, f)
        # throw a warning
        warnings.warn("Ran in no-run mode")
        return None, None
    # create number of reactors
    reactors = []
    integration_success = []
    for i in range(n_pz):
        # reset fuel object
        fuel.TP = T_i, P_i
        fuel.set_equivalence_ratio(phis[i], X_fuel, X_air, basis="mole")
        # create reservoirs
        inlet = ct.Reservoir(fuel)
        outlet = ct.Reservoir(fuel)
        # create reactor
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fuel.equilibrate("HP")
        reactors.append(ct.IdealGasConstPressureMoleReactor(fuel))
        reactors[i].volume = split_volumes[i]
        tres = split_volumes[i] / reactors[i].thermo.density / mfs[i]
        # connect
        mfc = ct.MassFlowController(inlet, reactors[i], mdot=lambda t: mfs[i])
        ct.PressureController(reactors[i], outlet, primary=mfc)
        # create a network and integrate the reactor to residence time
        # it is faster and more stable to integrate them all separately than together
        net = new_network([reactors[i]])
        net.max_time_step = 0.1
        startState = reactors[i].thermo.TPX
        # setup some restarts for failed cases
        success = False
        for j in range(3):
            try:
                net.advance(100 * tres)
                success = True
                break
            except Exception as e:
                print(e)
                print(f"Failed PZ integration, reducing timestep from {net.max_time_step}...")
                reactors[i].thermo.TPX = startState
                reactors[i].syncState()
                reactors[i].volume = split_volumes[i]
                net.initial_time = 0 # reset time to 0
                net.max_time_step = net.max_time_step / 5
                net.initialize()
        # add all successful integration attempts
        integration_success.append(success)
    # go through any unsuccessful attempts
    for p, sv in enumerate(integration_success):
        if not sv:
            print(f"Reactor case {p} failed...")
            # turn off preconditioning in unsuccessful cases
            try:
                net = new_network([reactors[p]], precon=False)
                net.max_time_step = 1e-4
                tres = split_volumes[p] / reactors[p].thermo.density / mfs[p]
                net.advance(100 * tres)
                integration_success[p] = True
            except Exception as e:
                raise Exception("Complete primary zone integration failure.")
    thermo_states["tres"] = float(tres)
    # EI closure
    def find_EI(r, i, sp):
        return r.Y[r.component_index(sp)] * mfs[i] * 1000 / (mdot_fuel * mass_flow_fractions[i])
    # record EI of CO and of NOx
    if "CO" in fuel.species_names:
        ei_co = []
        for i, r in enumerate(reactors):
            ei_co.append(float(find_EI(r, i, "CO")))
        thermo_states["EI_CO_pz"] = ei_co
    if "NO" and "NO2" in fuel.species_names:
        ei_nox = []
        for i, r in enumerate(reactors):
            cei = find_EI(r, i, "NO")
            cei += find_EI(r, i, "NO2")
            ei_nox.append(float(cei))
        thermo_states["EI_NOx_pz"] = ei_nox
    # temperature distribution
    thermo_states["temperature_distribution"] = [float(r.T) for r in reactors]
    # mix together at constant enthalpy and pressure
    mds = numpy.zeros(fuel.n_species)
    enthalpy = 0
    for i, r in enumerate(reactors):
        mds += r.Y * mfs[i]
        enthalpy += r.thermo.HP[0] * mfs[i]
        pressure = r.thermo.HP[1]
    mds = mds / mdot
    enthalpy /= mdot
    fuel.HPY = (enthalpy, pressure, mds)
    thermo_states["sz_init_phi"] = float(fuel.equivalence_ratio(X_fuel, X_air))
    print_state_TP(fuel.TP[0], fuel.TP[1], f"2pz:{thrust_level:0.2f}")
    thermo_states.get("T", []).append(f"{fuel.TP[0]:0.2f}")
    thermo_states.get("P", []).append(f"{fuel.TP[1]:0.2f}")
    # Instead of modeling as a PFR, converting dilution as a function of z
    # to dilution as a function of time
    # reactors
    fast_mixing = DilutionReactor(fuel, mdot=mdot/2, beta_da=beta_da, beta_mixing=beta_sa_fm, mixing_scale=l_sa_fm)
    slow_mixing = DilutionReactor(fuel, mdot=mdot/2, beta_da=beta_da, beta_mixing=beta_sa_sm, mixing_scale=l_sa_sm)
    # create dilution reservoir
    fuel.TPX = T_i, P_i, X_air
    # set needed dilution values
    for r in [fast_mixing, slow_mixing]:
        r.mass_fractions_air = fuel.Y
        r.partial_enthalpy_air = fuel.partial_molar_enthalpies * fuel.molecular_weights
        r.enthalpy_mass_air = fuel.enthalpy_mass
    # create network
    net = new_network([slow_mixing, fast_mixing], precon=False)
    net.initialize()
    # solution array for combustor
    if kwargs.get("output_mixing", False):
        fast_mixing_states = ct.SolutionArray(fast_mixing.thermo, extra=["z", "mdot", "t", "phi"])
        slow_mixing_states = ct.SolutionArray(slow_mixing.thermo, extra=["z", "mdot", "t", "phi"])
        while fast_mixing.zloc < L_sz or slow_mixing.zloc < L_sz:
            slow_mixing_states.append(slow_mixing.thermo.state, z=slow_mixing.zloc, t=0.0, mdot=slow_mixing.mass_flow_rate, phi=slow_mixing.thermo.equivalence_ratio(X_fuel, X_air))
            fast_mixing_states.append(fast_mixing.thermo.state, z=fast_mixing.zloc, t=0.0, mdot=fast_mixing.mass_flow_rate, phi=fast_mixing.thermo.equivalence_ratio(X_fuel, X_air))
            # take 10 steps before every record
            for i in range(50):
                net.step()
        # write last state
        slow_mixing_states.append(slow_mixing.thermo.state, z=slow_mixing.zloc, t=0.0, mdot=slow_mixing.mass_flow_rate, phi=slow_mixing.thermo.equivalence_ratio(X_fuel, X_air))
        fast_mixing_states.append(fast_mixing.thermo.state, z=fast_mixing.zloc, t=0.0, mdot=fast_mixing.mass_flow_rate, phi=fast_mixing.thermo.equivalence_ratio(X_fuel, X_air))
        # write out state and reformat
        fm5 = os.path.join(outdir, f"fast-mixing-states{name_id}.hdf5")
        fast_mixing_states.save(fm5, overwrite=True, name="thermo")
        reformat_hdf5(fm5)
        # write out state and reformat
        sm5 = os.path.join(outdir, f"slow-mixing-states{name_id}.hdf5")
        slow_mixing_states.save(sm5, overwrite=True, name="thermo")
        reformat_hdf5(sm5)
    else:
        while fast_mixing.zloc < L_sz or slow_mixing.zloc < L_sz:
            for i in range(10):
                net.step()
    # mix together at constant enthalpy and pressure
    zones = [slow_mixing, fast_mixing]
    mds = numpy.zeros(fuel.n_species)
    enthalpy = 0
    total_mdot = 0
    for i, r in enumerate(zones):
        total_mdot += r.mass_flow_rate
        mds += r.Y * r.mass_flow_rate
        enthalpy += r.thermo.HP[0] * r.mass_flow_rate
        pressure = r.thermo.HP[1]
    mds = mds / total_mdot
    enthalpy /= total_mdot
    fuel.HPY = (enthalpy, pressure, mds)
    print_state_TP(fuel.TP[0], fuel.TP[1], f"2sz:{thrust_level:0.2f}")
    thermo_states.get("T", []).append(f"{fuel.TP[0]:0.2f}")
    thermo_states.get("P", []).append(f"{fuel.TP[1]:0.2f}")
    thermo_states["mdot_out"] = float(slow_mixing.mass_flow_rate + fast_mixing.mass_flow_rate)
    # make and return combustor reactor
    combustor = ct.IdealGasConstPressureMoleReactor(fuel)
    return combustor, thermo_states


def surrogate_selection(choice=0):
    if (choice == 0):
        xfuel = "N-C10H22:0.4267, I-C8H18:0.3302, C7H8:0.2431"
    elif (choice == 1):
        xfuel = "N-C12H26:0.303, C7H14:0.485, ebz:0.212"
    elif (choice == 2):
        xfuel = "N-C12H26:0.384, C7H14:0.234, N-C16H34:0.148, C7H8:0.234"
    elif (choice == 3):
        xfuel = "N-C12H26:0.290, N-C16H34:0.142, C7H8:0.249, DECALIN:0.319"
    elif (choice == 4):
        xfuel = "N-C12H26:0.371, I-C8H18:0.020, N-C16H34:0.206, C7H8:0.259, DECALIN:0.145"
    elif (choice == 5):
        xfuel = "N-C12H26:0.404, I-C8H18:0.295, ebz:0.073, pbz:0.228"
    # napthalene calc
    items = {}
    sumt = 0
    for fc in xfuel.split(", "):
        n,v = fc.split(":")
        items[n] = float(v)
        sumt += float(v)
    if not numpy.isclose((1-sumt), 0) and 1-sumt > 0:
        items["C10H8"] = 1-sumt
    return items


@click.command()
@click.option('--equiv_ratio', default=1.0, help='Equivalence ratio for the fuel')
@click.option("--farnesane", default=0.1, help="Farnesane blend percentage")
@click.option("--outdir", default=DATA_DIR, help="Output directory for data")
@click.option("--thrust", default=1.0, help="Percentage of thrust to use")
@click.option("--stime", default=0.0, help="Time of day in hours when the simulation is started")
@click.option("--fmodel", default=os.path.join(DATA_DIR, "combustor-nox.yaml"), help="Fuel model used")
@click.option("--amodel", default=os.path.join(DATA_DIR, "atmosphere.yaml"), help="Atmospheric model used")
@click.option("--restdir", default=None, help="Restart directory for atmospheric only simulation")
@click.option('--precon_off', default=True, is_flag=True, help='Turn off preconditioner')
@click.option('--npz', default=21, help='Number of reactors used in simulation')
@click.option('--hightol', default=True, is_flag=True, help='Turn off preconditioner')
@click.option('--nox', default=0.0, help='Atmospheric NOX injection percentage, make value -1 for equilibrium estimate.')
@click.option('--h2o', default=0.0, help='Set air water content.')
@click.option('--name', default="", help='Append content to the name of files.')
@click.option('--pthermo', default=False, is_flag=True, help='Use PyCycle thermo data')
@click.option('--complex_ambient', default=False, is_flag=True, help='Use complex ambient estimation.')
@click.option('--entf', default=1.0, help='Entrainment rate scalar')
@click.option('--emodel', default="5B", help='Engine model, 7B or 5B')
@click.option('--total_phi', default=False, is_flag=True, help='Makes the specified equivalence ratio across the entire system, not the primary zone.')
@click.option('--surrogate', default=0, help='Number for surrogate fuel used.')
@click.option('--norun', default=False, is_flag=True, help='A testing flag to prevent full run.')
def run_combustor_atm_sim(equiv_ratio, farnesane, outdir, thrust, stime, fmodel, amodel, restdir, precon_off, npz, hightol, nox, h2o, name, pthermo, complex_ambient, entf, emodel, total_phi, surrogate, norun):
    global PRECONDITIONED
    PRECONDITIONED = precon_off
    global HIGH_TOLERANCE
    HIGH_TOLERANCE = hightol
    global ENGINE_MODEL
    ENGINE_MODEL = emodel
    # set  a new entrainment rate scalar
    global ENTRAINMENT_FACTOR
    ENTRAINMENT_FACTOR = entf
    PlumeReactor.entrainment_scalar = entf
    # run combustor sim
    combustor_atm_sim(equiv_ratio, farnesane, outdir, stime=stime, fmodel=fmodel, amodel=amodel, thrust_level=thrust, restart_dir=restdir, npz=npz, nox=nox, h2o=h2o, name=name, pthermo=pthermo, complex_ambient=complex_ambient, total_phi=total_phi, surrogate=surrogate, norun=norun)

def get_complex_ambient_by_model(solobj):
    with open(os.path.join(DATA_DIR, "complex-ambient.yaml"), "r") as caf:
        ca_ppm = yaml.load(caf, Loader=yaml.SafeLoader)
    total = 0
    all_names = solobj.species_names
    ca_ppm = dict(filter(lambda x: x[0] in all_names, ca_ppm.items()))
    for k, v in ca_ppm.items():
        ca_ppm[k] = float(v) / 1e6
        total += ca_ppm[k]
    ca_ppm["N2"] = 1 - total
    X_air = [f"{k}:{v}" for k, v in ca_ppm.items()]
    X_air = ", ".join(X_air)
    return X_air

def test_complex_ambient():
    # create path to models
    fuel_model = os.path.join(DATA_DIR, "combustor.yaml")
    # creation of fuel thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, name="combustor", transport_model=None)
    X_air = get_complex_ambient_by_model(fuel)

def combustor_atm_sim(equiv_ratio, farnesane, outdir, stime=0.0, fmodel=None, amodel=None, thrust_level=1.0, restart_dir=None, npz=21, nox=0.0, h2o=0.0, name="", pthermo=False, complex_ambient=False, entf=0, total_phi=False, surrogate=0, norun=False):
    # parameters used in both cases
    if name:
        name = f"-{name}"
    states_file = os.path.join(outdir, f"thermo-states-{equiv_ratio:.1f}-{farnesane:.2f}{name}.yaml")
    thermo_states = {"farnesane": f"{farnesane:.2f}", "equivalence_ratio": f"{equiv_ratio:.1f}", "engine_model":ENGINE_MODEL, "entrainment_rate_scalar":ENTRAINMENT_FACTOR}
    # attempt to make directory
    try:
        os.mkdir(outdir)
    except:
        pass
    # create path to models
    fuel_model = os.path.join(DATA_DIR, "combustor.yaml") if fmodel is None else fmodel
    atms_model = os.path.join(DATA_DIR, "atmosphere.yaml") if amodel is None else amodel
    # creation of fuel thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, name="combustor", transport_model=None)
        atms = PlumeSolution(atms_model, name="atmosphere", transport_model=None)
    # append appropriate directories
    sys.path.append(DATA_DIR)
    # create initial fuel setting
    X_fuel = surrogate_selection(surrogate)
    # adjust fuels for farnesane blend
    X_fuel = {k:v * (1-farnesane) for k,v in X_fuel.items()}
    # add farnesane blend
    Xfarne = 1 - numpy.sum([v for k,v in X_fuel.items()])
    if Xfarne > 0:
        X_fuel.update({"iC15H32": Xfarne})
    X_fuel = ", ".join([f"{k}:{v:.6f}" for k, v in X_fuel.items()])
    if complex_ambient:
        print("Using complex ambient configuration for air...")
        X_air = get_complex_ambient_by_model(fuel)
    else:
        X_air = f"H2O: {h2o:.2f}, O2:0.2095, N2:0.7808"
    thermo_states["X_fuel"] = X_fuel
    thermo_states["X_air"] =  X_air
    # design parameters
    p_atm = get_prop_by_thrust_level(thrust_level, "Pt0[kPa]") * 1000
    T_atm = get_prop_by_thrust_level(thrust_level, "Tt0[K]")
    print_state_TP(T_atm, p_atm, 0)
    thermo_states["T"] = [f"{T_atm:0.2f}"]
    thermo_states["P"] = [f"{p_atm:0.2f}"]
    # set initial conditions
    fuel.TPX = T_atm, p_atm, X_air
    atms.TPX = T_atm, p_atm, X_air
    # setup air properties needed later
    air_state = atms.TPX
    air_enthalpy = atms.partial_molar_enthalpies
    air_cp = atms.partial_molar_cp
    # restart from existing combustor state to perturb atms only
    if restart_dir is None:
        # compression is achieved in data file during multizone-combustor run
        combustor, ____ = multizone_combustor(fuel, thrust_level, equiv_ratio, X_fuel, X_air, thermo_states=thermo_states, outdir=outdir, name_id=f"-{equiv_ratio:.1f}-{farnesane:.2f}{name}", n_pz=npz, use_total_equiv=total_phi, norun=norun)
        if norun:
            return
        # continue
        T2, p2 = combustor.thermo.TP
        # get total moles for volume calculation
        Ntotal = numpy.sum(combustor.get_state()[1:fuel.n_species + 1])
        # use nozzle end conditions to back calculate state 3 and calculate T4
        A3 = numpy.pi * (get_prop_by_thrust_level(thrust_level, "D[in]") * 0.0254) ** 2 / 4
        A4 = numpy.pi * (get_prop_by_thrust_level(thrust_level, "ThroatDia[in]") * 0.0254) ** 2 / 4
        M4 = get_prop_by_thrust_level(thrust_level, "SecMach")
        p4 = p_atm
        # iteratively solve for pressure at 3 from choked nozzle and area ratio
        # to get temperature at 3 and then at 4
        M3 = 0.1 # initial guess
        gamma_old = 0
        gamma = fuel.cp_mole / (fuel.cp_mole - ct.gas_constant)
        while numpy.abs(gamma - gamma_old) > 1e-8:
            gmone = gamma - 1
            def machthree(x):
                numer = 1 + gmone / 2 * M4 * M4
                denom = 1 + gmone / 2 * x[0] * x[0]
                mach = A4 / A3 * M4 / ((numer / denom) ** ((gamma + 1) / (2 * gamma - 2)))
                return x[0] - mach
            M3 = fsolve(machthree, [M3])[0]
            # solve for pressure 3
            p3 = p4 / (((1 + gmone / 2 * M3 * M3) / (1 + gmone / 2 * M4 * M4)) ** (gamma / gmone))
            T3 = T2 * ((p3 / p2) ** (ct.gas_constant / fuel.cp_mole))
            # set new state and recalculate gamma
            fuel.TPY = T3, p3, combustor.Y
            gamma_old = gamma
            gamma = fuel.cp_mole / (fuel.cp_mole - ct.gas_constant)
        # print state at 3
        v3 = Ntotal * ct.gas_constant * T3 / p3
        print_state_TP(T3, p3, 3)
        print(f"M3: {M3}")
        T4 = T3 * (1 + gmone / 2 * (M3 ** 2)) / (1 + gmone / 2 * (M4 ** 2))
        v4 = Ntotal * ct.gas_constant * T4 / p4
        print_state_TP(T4, p4, 4)
        # fill thermo states
        thermo_states["T"].append(f"{T3:0.2f}")
        thermo_states["P"].append(f"{p3:0.2f}")
        thermo_states["T"].append(f"{T4:0.2f}")
        thermo_states["P"].append(f"{p4:0.2f}")
        thermo_states["M3"] = [f"{M3:0.2f}"]
        thermo_states["gamma"] = gamma
        # open model mapping
        mapfile = f'{os.path.basename(fmodel).split(".")[0]}_to_{os.path.basename(amodel).split(".")[0]}.yaml'
        mapfile = os.path.join(os.path.dirname(fmodel), mapfile)
        if not os.path.exists(mapfile):
            map_models_wrapped(fmodel, amodel)
        with open(mapfile, "r") as f:
            mapping = yaml.load(f, Loader=yaml.SafeLoader)
        Y_mapped = numpy.zeros(atms.n_species)
        cstate = combustor.Y
        for csp, asp in mapping.items():
            Y_mapped[atms.species_index(asp)] = cstate[fuel.species_index(csp)]
    else:
        print("Restarting from existing combustor data")
        # load existing thermo states
        sf = os.path.join(restart_dir, os.path.basename(states_file))
        with open(sf, "r") as f:
            ydata = yaml.load(f, Loader=yaml.SafeLoader)
            thermo_states.update(ydata)
        # loading solution array
        rf = os.path.join(restart_dir, f"short-term-states-{equiv_ratio:.1f}-{farnesane:.2f}{name}.hdf5")
        short_data = h5py.File(rf, "r")
        T4 = short_data["T"][0]
        A4 = numpy.pi * (0.62 ** 2) / 4 # meters
        M4 = get_prop_by_thrust_level(thrust_level, "SecMach")
        p4 = p_atm
        gamma = thermo_states.get("gamma", 1.4)
        Y_mapped = numpy.zeros(atms.n_species)
        for i in range(atms.n_species):
            Y_mapped[i] = short_data[f"Y_{atms.species_name(i)}"][0]
    # Use pycycle thermo conditions
    if pthermo:
        p4 = get_prop_by_thrust_level(thrust_level, "Pt5[kPa]") * 1000
        T4 = get_prop_by_thrust_level(thrust_level, "Tt5[K]")
    # relevant mach numbers
    thermo_states["M4"] = [f"{M4:0.2f}"]
    # inject nox or equilibrium estimate it
    if nox > 0:
        print(f"Injecting {nox} of NO2")
        Y_mapped[atms.species_index("NO2")] = nox
    elif nox == -1:
        nox_sp = ["NO", "NO2"]
        # reset fuel object
        P_i = get_prop_by_thrust_level(thrust_level, "Pt3[kPa]") * 1000
        T_i = get_prop_by_thrust_level(thrust_level, "Tt3[K]")
        fuel.TP = T_i, P_i
        fuel.set_equivalence_ratio(equiv_ratio, X_fuel, X_air, basis="mole")
        # create reactor
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fuel.equilibrate("HP")
        for nsp in nox_sp:
            print(f"Equilibrium estimate for {nsp}: {fuel.Y[fuel.species_index(nsp)]}")
            Y_mapped[atms.species_index(nsp)] = fuel.Y[fuel.species_index(nsp)]
    # create atmosphere reactor
    atmosphere = PlumeReactor(atms, start_time=stime)
    # set thermo object of atmospheric object and sync state
    atmosphere.thermo.TPY = T4, p4, Y_mapped
    atmosphere.syncState()
    # set atmospheric quantities
    atmosphere.TX_air = [air_state[0],] + [x for x in air_state[2]]
    atmosphere.enthalpy_air = air_enthalpy
    atmosphere.cp_air = air_cp
    # setup mass flow from reservoir to atmosphere
    gc = 1 # 1 in SI system
    mdot = atmosphere.thermo.density * A4 * M4 * numpy.sqrt(gamma * gc * 287 * T4)
    thermo_states["mdot_exhaust"] = float(mdot)
    # setup reactor network
    net = new_network([atmosphere], precon=True)
    # setup solution array
    short_states = ct.SolutionArray(atms, extra=["time", "moles", "mass", "volume"])
    failures = 0
    moles = numpy.sum(atmosphere.get_state()[1:])
    short_states.append(atmosphere.thermo.state, time=net.time, moles=moles, mass=atmosphere.mass, volume=atmosphere.volume)
    long_states = ct.SolutionArray(atms, extra=["time", "moles", "mass", "volume"])
    # setup times to safely integrate too and store data for because storing mass model
    # species data is expensive
    if not numpy.isclose(atmosphere.entrainment_scalar, 1.0):
        times = numpy.logspace(numpy.log10(1e-6), numpy.log10(500), 500)
    else:
        times = numpy.logspace(numpy.log10(1e-6), numpy.log10(10), 100)
    failed = False
    for t in times: #for t in steps:
        failures = 0
        t_target = t
        net.max_time_step = 0
        while net.time < t:
            try:
                print(f"Integrating atmosphere time, temp: {t_target, atmosphere.T}...")
                former_time = net.time
                net.advance(t_target)
                # reset targeted time and time step
                t_target = t
                failures = 0
                net.max_time_step = 0
                # append the state so last state is always available
                moles = numpy.sum(atmosphere.get_state()[1:atms.n_species + 1])
                short_states.append(atmosphere.thermo.state, time=net.time, mass=atmosphere.mass, moles=moles, volume=atmosphere.volume)
                long_states.append(atmosphere.thermo.state, time=net.time, mass=atmosphere.mass, moles=moles, volume=atmosphere.volume)
                # break if temperatures are close
                if check_entrainment_conditions(atmosphere):
                    print(f"Temperature equilibrium tolerance met, breaking...")
                    break
            except Exception as e:
                print(f"Failure {failures}: reducing max time step")
                # reset reactor and network
                atmosphere.thermo.TDX = short_states[-1].TDX
                atmosphere.syncState()
                atmosphere.volume = short_states[-1].volume
                net.initial_time = former_time
                if net.max_time_step == 0:
                    net.max_time_step = 0.1
                else:
                    net.max_time_step = net.max_time_step / 10
                net.initialize()
                t_target = former_time + t_target * 0.01 # 1% of next target time
                failures += 1
                if failures >= 20:
                    failed = True
                    break
        if failed or check_entrainment_conditions(atmosphere):
            break
    print(f"Integrated atmosphere time, temp: {net.time, atmosphere.T}...")
    # append thermo state
    thermo_states["T"].append(f"{atmosphere.T:0.2f}")
    thermo_states["P"].append(f"{p4:0.2f}")
    print_state_TP(atmosphere.T, p4, 5)
    # longer time scale integration - turn off entrainment and energy
    atmosphere.entrainment = False
    atmosphere.energy_enabled = False
    atmosphere.no_change_species = ["H2O", "O2", "N2"]
    failures = 0
    # setup times to safely integrate too and store data for because storing mass model
    # species data is expensive
    ndays = 3
    times = numpy.linspace(net.time, ndays * 24 * 3600, 100)
    failed = False
    for t in times[1:]:
        failures = 0
        t_target = t
        net.max_time_step = 1000
        while net.time < t:
            try:
                print(f"Integrating atmosphere time (days), temp: {net.time / (3600 * 24), atmosphere.T}...")
                former_time = net.time
                net.advance(t_target)
                # reset targeted time and time step
                t_target = t
                failures = 0
                net.max_time_step = 1000
                # append the state so last state is always available
                moles = numpy.sum(atmosphere.get_state()[1:atms.n_species + 1])
                long_states.append(atmosphere.thermo.state, time=net.time, mass=atmosphere.mass, moles=moles, volume=atmosphere.volume)
            except Exception as e:
                print(f"Failure {failures}: reducing max time step")
                # reset reactor and network
                atmosphere.thermo.TDX = long_states[-1].TDX
                atmosphere.syncState()
                atmosphere.volume = long_states[-1].volume
                net.initial_time = former_time
                net.max_time_step = net.max_time_step / 10
                net.initialize()
                t_target = former_time + t_target * 0.01 # 1% of next target time
                failures += 1
                if failures >= 20:
                    failed = True
                    break
        if failed:
            break
    # write out thermo yaml
    with open(states_file, "w") as f:
        yaml.dump(thermo_states, f)
    # write out hdf5 data - short term
    try:
        ah5 = os.path.join(outdir, f"short-term-states-{equiv_ratio:.1f}-{farnesane:.2f}{name}.hdf5")
        short_states.save(ah5, overwrite=True, name="thermo", sub="data")
        reformat_hdf5(ah5)
    except Exception as e:
        os.remove(ah5)
        print(f"Failed to save {ah5}, writing to csv...")
        ah5 = os.path.join(outdir, f"short-term-states-{equiv_ratio:.1f}-{farnesane:.2f}{name}.csv")
        short_states.save(ah5, overwrite=True)
    # write out hdf5 data - longterm
    try:
        ah5 = os.path.join(outdir, f"long-term-states-{equiv_ratio:.1f}-{farnesane:.2f}{name}.hdf5")
        long_states.save(ah5, overwrite=True, name="thermo")
        reformat_hdf5(ah5)
    except Exception as e:
        os.remove(ah5)
        print(f"Failed to save {ah5}, writing to csv...")
        ah5 = os.path.join(outdir, f"long-term-states-{equiv_ratio:.1f}-{farnesane:.2f}{name}.csv")
        long_states.save(ah5, overwrite=True)

def check_entrainment_conditions(r):
    val = True
    for i, v in enumerate(r.TX_air[1:]):
        if v > 0:
            val = numpy.allclose(r.thermo.X[i], v, rtol=1e-1) and val
    return numpy.isclose(r.T, r.TX_air[0], rtol=0.025) and val

def curve_fit_thrust_data(xval=0, mdot=1.086, test=False):
    with open(os.path.join(DATA_DIR, "thrust-data.csv"), "r")  as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)
        data = list(map(lambda x: (float(x[0]), float(x[1])), reader))
    xd, yd = zip(*data) # thrust level, mfraction
    yd = numpy.array(yd) / numpy.amax(yd) * mdot
    def func(x, a, b):
        return a + b * numpy.log(x) + 0.103
    popt, pcov = curve_fit(func, xd, yd)
    if test:
        x = numpy.linspace(0.01, 1.0, 10)
        y = [func(xi, popt[0], popt[1]) for xi in x]
        plt.plot(x, y)
        plt.show()
        print(x, y)
    return func(xval, popt[0], popt[1])


if __name__ == "__main__":
    run_combustor_atm_sim()
