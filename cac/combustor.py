import os
import sys
import h5py
import click
import numpy
import warnings
import cantera as ct
import yaml
from scipy.optimize import fsolve
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
            hf.create_dataset(f"Y_{c}", data=Y_arr[:, i])
        # write the rest of the data out
        other_keys = list(arr.keys())
        other_keys.remove("Y")
        for n in other_keys:
            data = hf["thermo"]["data"][n][()]
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

def braggs_combustor(fuel, thermo_states={"T":[], "P":[]}):
    # combustor geometrical parameters
    plenum_radius = 0.035
    plenum_length = 0.1
    chamber_radius = 0.045
    chamber_length = 0.3
    # residence time in combustor
    residence_time = 1.0  # starting residence time
    # initial state
    initial_state = list(fuel.TPX)
    # inlet fuel tank
    fuel_tank = ct.Reservoir(fuel)
    # creating braggs combustor with a WSR and PFR
    fuel.equilibrate("HP")
    wsr = ct.IdealGasConstPressureMoleReactor(fuel)
    wsr.volume = numpy.pi * plenum_radius * plenum_radius * plenum_length
    # outlet exhaust reservoir
    combustion_chamber = ct.Reservoir(fuel)
    def mflow(t):
        return wsr.mass / residence_time
    # connect to reservoirs and get steady state operation
    ft_to_wsr = ct.MassFlowController(fuel_tank, wsr, mdot=mflow)
    wsr_to_cc = ct.PressureController(fuel_tank, combustion_chamber, primary=ft_to_wsr, K=0.01)
    net = new_network([wsr])
    # advance to steady state
    states = ct.SolutionArray(fuel, extra=['tres', "time", 'mass'])
    last_good_rt = residence_time
    failures = 0
    T_jump = 0
    net.max_time_step = 1e-6
    states.append(wsr.thermo.state, tres=residence_time, time=net.time, mass=wsr.mass)
    # loop to get starting temp
    while T_jump < 200 and wsr.T > 800:
        net.initial_time = 0.0  # reset the integrator
        try:
            T_f = wsr.T
            net.advance_to_steady_state(residual_threshold=1e-4)
            print('tres = {:.2e}, T = {:.1f}'.format(residence_time, wsr.T))
            states.append(wsr.thermo.state, tres=residence_time, time=net.time, mass=wsr.mass)
            if net.max_time_step < 1e-6:
                net.max_time_step = 1e-6
            residence_time *= 0.9  # decrease the residence time for the next iteration
            failures = 0
        except Exception as e:
            failures += 1 # increase number of failures
            net.max_time_step = net.max_time_step / 10 # decrease
            # reset state
            wsr.thermo.TPX = states[-1].TPX
            wsr.syncState()
            residence_time = states[-1].tres * 0.95 # reduce the residence time by smaller amount
            print('Failure: tres = {:.2e}; T = {:.1f}'.format(residence_time, wsr.T))
            if (failures >= 5): # 5 consecutive failures
                break
        T_jump = T_f - wsr.T
    # get a state a couple before the last
    initial_state[0] = states[-2].TP[0]
    fuel.TPX = initial_state
    residence = states[-2].tres
    thermo_states["T"].append(f"{fuel.TP[0]:0.2f}")
    thermo_states["P"].append(f"{fuel.TP[1]:0.2f}")
    print_state_TP(fuel.TP[0], fuel.TP[1], "2i")
    # create the combustor that acts as a plug flow reactor
    combustor = ct.IdealGasConstPressureMoleReactor(fuel)
    combustor.volume = numpy.pi * chamber_radius * chamber_radius * chamber_length
    # solution array for combustor
    combustor_states = ct.SolutionArray(combustor.thermo, extra=["U", "z", "t"])
    net = new_network([combustor])
    net.atol = 1e-24
    # WSR portion
    steady_state_advance(net, combustor, combustor_states)
    wsr_time = net.time
    wsr_len = len(combustor_states)
    # plug flow reactor states
    n_steps_pfr = 1000
    t_total = residence_time
    dt = t_total / n_steps_pfr
    mdot0 = combustor.mass / residence_time
    combustor_area = numpy.pi * chamber_radius * chamber_radius
    # define time, space, and other information vectors
    t1 = (numpy.arange(n_steps_pfr) + 1) * dt + net.time
    z1 = numpy.zeros_like(t1)
    combustor_states.append(combustor.thermo.state, U=mdot0 / combustor_area / combustor.thermo.density, z=0.0, t=net.time)
    for n1, t_i in enumerate(t1):
        # perform time integration
        net.advance(t_i)
        # compute velocity and transform into space
        U = mdot0 / combustor_area / combustor.thermo.density
        z1[n1] = z1[n1 - 1] + U * dt
        combustor_states.append(combustor.thermo.state, U=U, z=z1[n1], t=t_i)
    return combustor, combustor_states, wsr, states, (wsr_time, wsr_len)

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

if __name__ == "__main__":
    run_combustor_atm_sim()
