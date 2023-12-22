import os
import sys
import csv
import h5py
import yaml
import copy
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
from cac.constants import DATA_DIR
from cac.reactors import PlumeSolution, PlumeReactor

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


def new_network(r):
    # setup reactor network
    net = ct.ReactorNet(r)
    net.derivative_settings = {"skip-falloff": True,
                            "skip-third-bodies": True,
                            "skip-coverage-dependence": True,
                            "skip-electrochemistry": True,
                            "skip-flow-devices": True,
                            "skip-walls": True}
    net.max_steps = 100000
    # setup preconditioner
    net.preconditioner = ct.AdaptivePreconditioner()
    net.initialize()
    return net


@click.command()
@click.option('--equiv_ratio', default=1.0, help='Equivalence ratio for the fuel')
@click.option('--nsteps', default=100, help='The number of steps to take between records')
@click.option("--farnesane", default=0.1, help="Farnesane blend percentage")
@click.option("--sulfur", default=0.0001, help="Sulfur amount in mole fraction")
@click.option("--outdir", default=DATA_DIR, help="Output directory for data")
def run_combustor_atm_sim(equiv_ratio, nsteps, farnesane, sulfur, outdir):
    combustor_atm_sim(equiv_ratio, nsteps, farnesane, sulfur, outdir)


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


def get_mass_flow_thrust_data(thrust_level):
    cycle_data = pd.read_csv(os.path.join(DATA_DIR, "verification", 'CFM56-7B-737.csv'), delimiter=",", engine="python")
    ma = cycle_data['W3[kg/s]'].values
    mf = cycle_data['Wf[kg/s]'].values
    thrust = cycle_data['NetThrust[kN]']
    thrust_percents = thrust / numpy.amax(thrust)
    zipped = list(zip(thrust_percents, thrust, mf, ma))
    thrust_percents, thrust, mf, ma = zip(*sorted(zipped, key=lambda x: x[0]))
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


def multizone_combustor(fuel, thrust_level, equiv_pz, X_fuel, X_air, **kwargs):
    # default parameters
    thermo_states = kwargs.get("thermo_states", {"T":[], "P":[]})
    n_pz = kwargs.get("n_pz", 9)
    n_pz += 1 if n_pz % 2 == 0  else 0 # ensure always odd
    mixing_param = kwargs.get("mixing_param", 0.39) # mixing parameter
    volume = kwargs.get("volume", 0.0023) # m^3
    BPR = kwargs.get("BPR", 5.2)
    m_dot_fuel, m_dot_air, m_dot_fuel_takeoff, m_dot_air_takeoff = get_mass_flow_thrust_data(thrust_level)
    assert thrust_level <= 1.0
    equiv_sec_zone = kwargs.get("equiv_sz", 0.7)
    outdir = kwargs.get("outdir", "./")
    name_id = kwargs.get("name_id", "")
    # fixed params from study
    A_sz = 0.15 # m^2
    L_sz = 0.075 # m
    l_das = 0.95 # start
    l_dae = 1.0 # end
    f_sm = 0.5 # fraction of mass flow in slow mixing
    f_fm = 1 - f_sm
    l_sa_sm = 0.55
    l_sa_fm = 0.055
    T_i, P_i = fuel.TP
    # calculate total mass flow rate
    akeys = [fuel.species_index(x.split(":")[0].strip()) for x in X_air.split(",")]
    fkeys = [fuel.species_index(x.split(":")[0].strip()) for x in X_fuel.split(",")]
    # stoichiometric values
    fuel.set_equivalence_ratio(1.0, X_fuel, X_air)
    m_dot_fuel_st = m_dot_fuel_takeoff
    m_dot_st = numpy.sum(m_dot_fuel_st / fuel.Y[fkeys])
    m_dot_air_st= numpy.sum(m_dot_st * fuel.Y[akeys])
    FARst = m_dot_fuel_st / m_dot_air_st
    # fraction of total air
    f_air_pz = m_dot_fuel_takeoff / (m_dot_air_takeoff * equiv_pz * FARst)
    f_air_sa = m_dot_fuel_takeoff / (equiv_sec_zone * FARst * m_dot_air_takeoff)

    # update mass flow rates
    m_dot_air = m_dot_air * f_air_pz
    m_dot = m_dot_air + m_dot_fuel
    # now determine mean equiv ratio with takeoff air and mdot fuel
    equiv_mean = m_dot_fuel / (m_dot_air * FARst)
    # calculated parameters
    f_air_da = 1 - f_air_sa - f_air_pz
    beta_sa_sm = f_air_sa * f_sm * m_dot_air / (l_sa_sm * L_sz)
    beta_sa_fm = f_air_sa * f_fm * m_dot_air / (l_sa_fm * L_sz)
    beta_da = f_air_da * m_dot_air / (l_dae - l_das) / L_sz
    print(thrust_level, beta_sa_fm, beta_sa_sm, beta_da)
    # calculate phi values
    lower_bound = (0 - equiv_mean) / (mixing_param * equiv_mean)
    norm_dist = stats.truncnorm(lower_bound, 5, loc=equiv_mean, scale=mixing_param * equiv_mean)
    phis = numpy.linspace(norm_dist.ppf(0.01), norm_dist.ppf(0.99), n_pz)
    thermo_states["phi_distribution"] = [float(p) for p in phis]
    # calculate mass flows
    mass_flow_fractions = (1 / (mixing_param * numpy.sqrt(2*numpy.pi))) * numpy.exp(- (phis - equiv_mean) ** 2 / (2 * mixing_param ** 2 )) * (phis[1] - phis[0])
    mfs = mass_flow_fractions * m_dot
    # add to thermo states
    thermo_states["mass_flow_fraction_distribution"] = [float(p) for p in mass_flow_fractions]
    thermo_states["mass_flow_distribution"] = [float(p) for p in mfs]
    split_volume = volume / n_pz # split over zones
    # create number of reactors
    reactors = []
    initial_mass_fuel = 0
    for i in range(n_pz):
        # reset fuel object
        fuel.TP = T_i, P_i
        fuel.set_equivalence_ratio(phis[i], X_fuel, X_air, basis="mole")
        # create reservoirs
        inlet = ct.Reservoir(fuel)
        outlet = ct.Reservoir(fuel)
        # create reactor
        reactors.append(ct.IdealGasConstPressureMoleReactor(fuel))
        reactors[i].volume = split_volume
        # get mass of initial fuel for verification purposes
        initial_mass_fuel += numpy.sum(reactors[i].Y[fkeys] * reactors[i].mass)
        reactors[i].thermo.equilibrate("HP")
        reactors[i].syncState()
        # connect
        mfc = ct.MassFlowController(inlet, reactors[i], mdot=mfs[i])
        ct.PressureController(reactors[i], outlet, primary=mfc)
    # create network and advance to steady state to simulate well stirred reactors
    net = new_network(reactors)
    # net.max_time_step = 1e-7
    net.advance_to_steady_state(residual_threshold=1e-2)
    # temperature distribution
    thermo_states["temperature_distribution"] = [r.T for r in reactors]
    thermo_states["initial_fuel_mass"] = float(initial_mass_fuel)
    # mix together at constant enthalpy and pressure
    quantities = [ct.Quantity(r.thermo, constant="HP") for r in reactors]
    Q = quantities[0]
    for q in quantities[1:]:
        Q += q
    fuel.TPX = Q.TPX
    # print_state_TP(fuel.TP[0], fuel.TP[1], "2pz")
    thermo_states["T"].append(fuel.TP[0])
    thermo_states["P"].append(fuel.TP[1])
    # PFR representing the rest of the combustor
    # reactors
    fast_mixing = ct.IdealGasConstPressureMoleReactor(fuel)
    fast_mixing.volume = volume / 2
    slow_mixing = ct.IdealGasConstPressureMoleReactor(fuel)
    slow_mixing.volume = volume / 2
    # adjust volume to preserve mass
    total_mass = numpy.sum([r.mass for r in reactors])
    adjust_volume_to_preserve_mass(fast_mixing, total_mass / 2)
    slow_mixing.volume = fast_mixing.volume
    assert numpy.isclose(fast_mixing.mass, total_mass / 2)
    assert numpy.isclose(slow_mixing.mass, total_mass / 2)
    # setup PFR
    n_steps = 1000
    u = m_dot / (fast_mixing.density * A_sz)
    t_total = L_sz / u
    dt = t_total / n_steps
    # dilution function for slow zone
    Lcs = 0
    dilution_on = True
    def slow_zone_dilution(t):
        if dilution_on:
            if Lcs >= 0 and Lcs <= l_sa_sm * L_sz:
                return beta_sa_sm * Lcs
            elif Lcs >= l_das * L_sz and Lcs <= l_dae * L_sz:
                return beta_da * Lcs
        return 0
    # dilution function for fast zone
    Lcf = 0
    def fast_zone_dilution(t):
        if dilution_on:
            if Lcf >= 0 and Lcf <= l_sa_fm * L_sz:
                # print("Lcf", Lcf)
                return beta_sa_fm * Lcf
            elif Lcf >= l_das * L_sz and Lcf <= l_dae * L_sz:
                return beta_da * Lcf
        return 0

    # create dilution reservoir
    fuel.TP = T_i, P_i
    fuel.set_equivalence_ratio(0, X_fuel, X_air)
    dilution_air = ct.Reservoir(fuel)
    mfc_sm = ct.MassFlowController(dilution_air, slow_mixing, mdot=slow_zone_dilution)
    mfc_fm = ct.MassFlowController(dilution_air, fast_mixing, mdot=fast_zone_dilution)
    # create network
    netf = new_network([fast_mixing,])
    # vector of locations
    zf = [0.0]
    uf = [u]
    # integrate over length
    ctr = 1
    # solution array for combustor
    fast_mixing_states = ct.SolutionArray(fast_mixing.thermo, extra=["U", "z", "t"])
    fast_mixing_states.append(fast_mixing.thermo.state, U=uf[0], z=zf[0], t=0.0)
    while Lcf <= L_sz:
        # perform time integration
        netf.advance(ctr * dt)
        # compute velocity and transform into space
        uf.append(m_dot / A_sz / fast_mixing.thermo.density)
        # locations
        zf.append(zf[ctr - 1] + uf[ctr] * dt)
        Lcf = zf[ctr]
        # update array
        fast_mixing_states.append(fast_mixing.thermo.state, U=uf[ctr], z=zf[ctr], t=ctr*dt)
        ctr += 1
    # write out state and reformat
    fm5 = os.path.join(outdir, f"fast-mixing-states{name_id}.hdf5")
    fast_mixing_states.save(fm5, overwrite=True, name="thermo")
    reformat_hdf5(fm5)
    # create network
    nets = new_network([slow_mixing,])
    # vector of locations
    zs = [0.0]
    us = [u]
    # integrate over length
    ctr = 1
    # solution array for combustor
    slow_mixing_states = ct.SolutionArray(slow_mixing.thermo, extra=["U", "z", "t"])
    slow_mixing_states.append(slow_mixing.thermo.state, U=uf[0], z=zf[0], t=0.0)
    while Lcs <= L_sz:
        # perform time integration
        nets.advance(ctr * dt)
        # compute velocity and transform into space
        us.append(m_dot / A_sz / slow_mixing.thermo.density)
        # locations
        zs.append(zs[ctr - 1] + us[ctr] * dt)
        Lcs = zs[ctr]
        # update array
        slow_mixing_states.append(slow_mixing.thermo.state, U=us[ctr], z=zs[ctr], t=ctr*dt)
        ctr += 1
    # write out state and reformat
    sm5 = os.path.join(outdir, f"slow-mixing-states{name_id}.hdf5")
    slow_mixing_states.save(sm5, overwrite=True, name="thermo")
    reformat_hdf5(sm5)
    # mix together at constant enthalpy and pressure
    zones = [slow_mixing, fast_mixing]
    quantities = [ct.Quantity(r.thermo, constant="HP") for r in zones]
    Q = quantities[0]
    for q in quantities[1:]:
        Q += q
    fuel.TPX = Q.TPX
    # print_state_TP(fuel.TP[0], fuel.TP[1], "2sz")
    thermo_states["T"].append(fuel.TP[0])
    thermo_states["P"].append(fuel.TP[1])
    # make and return combustor reactor
    combustor = ct.IdealGasConstPressureMoleReactor(fuel)
    total_mass = slow_mixing.mass + fast_mixing.mass
    combustor.volume = slow_mixing.volume + fast_mixing.volume
    adjust_volume_to_preserve_mass(combustor, total_mass)
    assert numpy.isclose(combustor.mass, total_mass)

    return combustor, thermo_states
    # return (0, 0)


def combustor_atm_sim(equiv_ratio, nsteps, farnesane, sulfur, outdir, fmodel="combustor-min.yaml"):
    thermo_states = {"farnesane": f"{farnesane:.2f}", "sulfur": f"{sulfur:.4f}", "equivalence_ratio": f"{equiv_ratio:.1f}"}
    # append appropriate directories
    sys.path.append(DATA_DIR)
    # create initial fuel setting
    X_fuel = {"N-C10H22":0.4267, "I-C8H18":0.3302, "C7H8":0.2431}
    # add H2S
    if sulfur > 0:
        X_fuel.update({"H2S":sulfur})
    # adjust fuels for farnesane blend
    X_fuel = {k:v * (1-farnesane-sulfur) for k,v in X_fuel.items()}
    # add farnesane blend
    Xfarne = 1 - numpy.sum([v for k,v in X_fuel.items()])
    if Xfarne > 0:
        X_fuel.update({"iC15H32": Xfarne})
    X_fuel = ", ".join([f"{k}:{v:.6f}" for k, v in X_fuel.items()])
    X_air = "H2O: 0.04, O2:0.2095, N2:0.7808"
    thermo_states["X_fuel"] = X_fuel
    thermo_states["X_air"] =  X_air
    # design parameters
    EPR = 25 # engine pressure ratio
    EGT = 800 # K, exhaust gas temperatures
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    print_state_TP(T_atm, p_atm, 0)
    thermo_states["T"] = [f"{T_atm:0.2f}"]
    thermo_states["P"] = [f"{p_atm:0.2f}"]
    # create path to fuel model
    fuel_model = os.path.join(DATA_DIR, fmodel)
    # creation of fuel thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, name="combustor", transport_model=None, basis="mole")
    fuel.TPX = T_atm, p_atm, X_air
    # isentropic compression
    p1 = p_atm * EPR
    T1 = T_atm * (p1 / p_atm) ** (ct.gas_constant / fuel.cp_mole)
    thermo_states["T"].append(f"{T1:0.2f}")
    thermo_states["P"].append(f"{p1:0.2f}")
    print_state_TP(T1, p1, 1)
    # set input fuel conditions
    fuel.TP = T1, p1
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, X_air, basis="mole")
    # run braggs combustor
    combustor, combustor_states, wsr, wsr_states. ss_data = braggs_combustor(fuel, thermo_states=thermo_states)
    # write data
    ch5 = os.path.join(outdir, "combustor", f"combustor-pfr-states-{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}.hdf5")
    combustor_states.save(ch5, overwrite=True, name="thermo")
    reformat_hdf5(ch5)
    # continue
    T2, p2 = combustor.thermo.TP
    print_state_TP(T2, p2, 2)
    thermo_states["T"].append(f"{T2:0.2f}")
    thermo_states["P"].append(f"{p2:0.2f}")
    # get total moles for volume calculation
    moles = combustor.get_state()[1:]
    clean = True
    if clean:
        max_mole_value = numpy.amax(moles)
        decimals = int(numpy.ceil(numpy.abs(numpy.log10(max_mole_value))))
        tolerance = float(f"1e-{decimals+16}")
        moles = [m if m >= tolerance else 0 for m in moles]
    Ntotal = numpy.sum(moles)
    # use nozzle end conditions to back calculate state 3 and calculate T4
    A3 = numpy.pi * (1.07 ** 2) / 4 # meters
    A4 = numpy.pi * (0.62 ** 2) / 4 # meters
    M4 = 1
    p4 = p_atm
    # iteratively solve for pressure at 3 from choked nozzle and area ratio
    # to get temperature at 3 and then at 4
    M3 = 0.3 # initial guess
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
    # create atmosphere model
    exhaust = ct.Reservoir(fuel)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        atms = PlumeSolution(fuel_model, name="atmosphere", transport_model=None, basis="mole")
    atms.TPX = T_atm, p_atm, X_air
    air_state = atms.TPX
    air_enthalpy = atms.partial_molar_enthalpies
    air_cp = atms.partial_molar_cp
    # create atmosphere reactor
    X_moles = moles / Ntotal
    atms.TPY = T4, p4, X_moles
    atmosphere = PlumeReactor(atms, energy="off")
    atmosphere.volume = v4
    atmosphere.TX_air = [air_state[0],] + [x for x in air_state[2]]
    atmosphere.enthalpy_air = air_enthalpy
    atmosphere.cp_air = air_cp
    # fill thermo states
    thermo_states["T"].append(f"{T3:0.2f}")
    thermo_states["P"].append(f"{p3:0.2f}")
    thermo_states["T"].append(f"{T4:0.2f}")
    thermo_states["P"].append(f"{p4:0.2f}")
    # setup mass flow from reservoir to atmosphere
    gc = 1 # 1 in SI system
    mdot = p4 * A4 * M4 * numpy.sqrt(gamma * gc / ct.gas_constant * T4)
    # setup reactor network
    net = new_network([atmosphere])
    net.max_time_step = 1e-6
    onsteps = nsteps
    # setup solution array
    atms_states = ct.SolutionArray(atms, extra=["time", "moles", "mass"])
    failures = 0
    moles = numpy.sum(atmosphere.get_state()[1:])
    atms_states.append(atmosphere.thermo.state, time=net.time, moles=moles, mass=atmosphere.mass)
    while atmosphere.T > T_atm + 10:
        print(f"Integrating atmosphere time, temp: {net.time, atmosphere.T}...")
        try:
            for i in range(nsteps):
                net.step()
                if atmosphere.T > 300:
                    break
            failures = 0
            if net.max_time_step < 1e-6:
                net.max_time_step = 1e-6
                nsteps = onsteps
        except Exception as e:
            failures += 1
            net.max_time_step = net.max_time_step / 10
            nsteps *= 10
            # reset state
            atmosphere.thermo.TPX = atms_states[-1].TPX
            atmosphere.syncState()
            net.initial_time = atms_states[-1].time
            if failures > 4:
                break
        # save state that makes it through
        moles = numpy.sum(atmosphere.get_state()[1:])
        atms_states.append(atmosphere.thermo.state, time=net.time, moles=moles, mass=atmosphere.mass)
    # append thermo state
    thermo_states["T"].append(f"{atmosphere.T:0.2f}")
    thermo_states["P"].append(f"{p4:0.2f}")
    print_state_TP(atmosphere.T, p4, 5)
    # write out thermo yaml
    states_file = os.path.join(outdir, "combustor", f"thermo-states-{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}.yaml")
    with open(states_file, "w") as f:
        yaml.dump(thermo_states, f)
    # write out hdf5 data
    ah5 = os.path.join(outdir, "combustor", f"atms-states-{equiv_ratio:.1f}-{farnesane:.2f}-{sulfur:.4f}.hdf5")
    atms_states.save(ah5, overwrite=True, name="thermo")
    reformat_hdf5(ah5)

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
