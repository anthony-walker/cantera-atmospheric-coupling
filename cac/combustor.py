import os
import sys
import click
import numpy
import warnings
import cantera as ct
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from cac.constants import DATA_DIR
from cac.reactors import PlumeSolution, PlumeReactor


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


def new_network(r, precon):
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
    if precon:
        net.preconditioner = ct.AdaptivePreconditioner()
    return net


@click.command()
@click.option('--equiv_ratio', default=1.0, help='Equivalence ratio for the fuel')
@click.option('--test', default=False, help='Simply just test integration without recording')
@click.option('--steadystate/--no-steadystate', default=False, help='Integrate directly to steady state')
@click.option('--precon/--no-precon', default=True, help='Add preconditioning to the simulation')
@click.option('--endtime', default=1.0, help='Run the simulation to the set end time')
@click.option('--nsteps', default=10, help='The number of steps to take between records')
@click.option("--farnesane", default=0.1, help="Farnesane blend percentage")
@click.option("--sulfur", default=0.001, help="Sulfur amount in mole fraction")
def run_combustor_atm_sim(equiv_ratio, test, steadystate, precon, endtime, nsteps, farnesane, sulfur):
    # append appropriate directories
    sys.path.append(DATA_DIR)
    # creation of combustor portion of simulation
    fuel_model = os.path.join(DATA_DIR, f"combustor-min.yaml")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, name="combustor", transport_model=None)
    # create initial fuel setting
    X_fuel = {"N-C10H22":0.4267, "I-C8H18":0.3302, "C7H8":0.2431}
    # add H2S
    if sulfur > 0:
        X_fuel.update({"H2S":sulfur})
    # adjust fuels for farnesane blend
    X_fuel = {k:v * (1-farnesane) for k,v in X_fuel.items()}
    # add farnesane blend
    Xfarne = 1 - numpy.sum([v for k,v in X_fuel.items()])
    if Xfarne > 0:
        X_fuel.update({"iC15H32": Xfarne})
    X_fuel = ", ".join([f"{k}:{v}" for k, v in X_fuel.items()])
    X_air = "H2O: 0.04, O2:0.2095, N2:0.7808"
    # combustor geometrical parameters
    plenum_radius = 0.035
    plenum_length = 0.1
    chamber_radius = 0.045
    chamber_length = 0.3
    # residence time in combustor
    residence_time = 0.1  # starting residence time
    # design parameters
    EPR = 25 # engine pressure ratio
    EGT = 800 # K, exhaust gas temperatures
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    fuel.TPX = T_atm, p_atm, X_air
    # isentropic compression
    p1 = p_atm * EPR
    T1 = T_atm * (p1 / p_atm) ** (ct.gas_constant / fuel.cp_mole)
    # set input fuel conditions
    fuel.TP = T1, p1
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, X_air)
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
    net = new_network([wsr], precon)
    # advance to steady state
    states = ct.SolutionArray(fuel, extra=['tres'])
    last_good_rt = residence_time
    failures = 0
    while wsr.T > 900:
        net.initial_time = 0.0  # reset the integrator
        try:
            net.advance_to_steady_state(residual_threshold=1e-2)
            print('tres = {:.2e}, T = {:.1f}'.format(residence_time, wsr.T))
            states.append(wsr.thermo.state, tres=residence_time)
            residence_time *= 0.9  # decrease the residence time for the next iteration
            failures = 0
        except:
            failures += 1
            last_known_rt = residence_time
            print('Failure: tres = {:.2e}; T = {:.1f}'.format(residence_time, wsr.T))
            if (failures > 10): # 10 consecutive failures
                break
            residence_time *= 0.99  # decrease the residence time for the next iteration
    # Plot results
    f, ax1 = plt.subplots(1, 1)
    ax1.plot(states.tres, states.heat_release_rate, '.-', color='C0')
    ax2 = ax1.twinx()
    ax2.plot(states.tres[:-1], states.T[:-1], '.-', color='C1')
    ax1.set_xlabel('residence time [s]')
    ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
    ax2.set_ylabel('temperature [K]', color='C1')
    f.tight_layout()
    plt.savefig(os.path.join(DATA_DIR, "combustor", f"heat-release-rate-{equiv_ratio}-{farnesane}-{sulfur}.pdf"))
    # get a state a couple before the last
    fuel.TDY = states[-2].TDY
    residence = states[-2].tres
    # create the combustor that acts as a plug flow reactor
    combustor = ct.IdealGasConstPressureMoleReactor(fuel)
    combustor.volume = numpy.pi * chamber_radius * chamber_radius * chamber_length
    net = new_network([combustor], precon)
    # plug flow reactor states
    n_steps_pfr = 1000
    t_total = residence_time
    dt = t_total / n_steps_pfr
    mdot0 = combustor.mass / residence_time
    combustor_area = numpy.pi * chamber_radius * chamber_radius
    # define time, space, and other information vectors
    t1 = (numpy.arange(n_steps_pfr) + 1) * dt
    z1 = numpy.zeros_like(t1)
    u1 = numpy.zeros_like(t1)
    combustor_states = ct.SolutionArray(combustor.thermo)
    for n1, t_i in enumerate(t1):
        # perform time integration
        net.advance(t_i)
        # compute velocity and transform into space
        u1[n1] = mdot0 / combustor_area / combustor.thermo.density
        z1[n1] = z1[n1 - 1] + u1[n1] * dt
        combustor_states.append(combustor.thermo.state)
    combustor_states.save(os.path.join(DATA_DIR, "combustor", f"combustor-pfr-states-{equiv_ratio}-{farnesane}-{sulfur}.csv"), overwrite=True)
    T2, p2 = combustor.thermo.TP
    # get total moles for volume calculation
    moles = combustor.get_state()[1:]
    max_mole_value = numpy.amax(moles)
    decimals = int(numpy.ceil(numpy.abs(numpy.log10(max_mole_value))))
    tolerance = float(f"1e-{decimals+16}")
    moles = [m if m >= tolerance else 0 for m in moles]
    Ntotal = numpy.sum(moles)
    # print state at 2
    print(T2, p2, combustor.volume)
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
    print(M3, T3, p3, v3)
    T4 = T3 * (1 + gmone / 2 * (M3 ** 2)) / (1 + gmone / 2 * (M4 ** 2))
    v4 = Ntotal * ct.gas_constant * T4 / p4
    print(T4, p4, v4)
    # create atmosphere model
    exhaust = ct.Reservoir(fuel)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        atms = PlumeSolution(fuel_model, name="atmosphere", transport_model=None)
    atms.TPX = T_atm, p_atm, X_air
    air_state = atms.TPX
    far_field = ct.Reservoir(atms)
    # create atmosphere reactor
    X_moles = moles / Ntotal
    atms.TPY = T4, p4, X_moles
    atmosphere = PlumeReactor(atms)
    atmosphere.volume = v4
    atmosphere.TX_air = [air_state[0],] + [x for x in air_state[2]]
    # create inlet reservoir for atmosphere
    exhaust = ct.Reservoir(atms)
    # setup mass flow from reservoir to atmosphere
    gc = 1 # 1 in SI system
    mdot = p4 * A4 * M4 * numpy.sqrt(gamma * gc / ct.gas_constant * T4)
    exhaust_mfc = ct.MassFlowController(exhaust, atmosphere, mdot=mdot)
    # setup reactor network
    net = new_network([atmosphere], precon)
    while atmosphere.T > T_atm:
        print(atmosphere.T)
        net.step()
    # atms_states = ct.SolutionArray(atms)
    # times = []
    # # adding the initial conditions
    # atms_states.append(atmosphere.thermo.state)
    # times.append(net.time)
    # # # loop for an hour of simulation time
    # # while net.time < endtime:
    # #     print_reactor_stats(net, atmosphere)
    # #     # try:
    # #     for i in range(nsteps):
    # #         net.step()
    # #         if net.time > endtime:
    # #             break
    # #     atms_states.append(atmosphere.thermo.state)
    # #     times.append(net.time)
    # #     # except Exception as e:
    # #     #     break
    # # print_reactor_stats(net, atmosphere)
    # # # investigations
    # # # C7H8, NA, SA
    # # # some plotting
    # # if not test:
    # #     # Plot results
    # #     f, ax1 = plt.subplots(1, 1)
    # #     ax1.plot(times, comb_states("CH4").Y*combustor.mass, '.-', color='b')
    # #     ax2 = ax1.twinx()
    # #     ax2.plot(times, atms_states("O3").Y*atmosphere.mass, '.-', color='g')
    # #     ax2.plot(times, atms_states("CO2").Y*atmosphere.mass, '.-', color='r')
    # #     # ax2.plot(times, atms_states("NO2").Y, '--', color='g')
    # #     # ax2.plot(times, atms_states("NO3").Y, '.-', color='g')
    # #     # ax2.plot(times, atms_states.T, '.-', color='C1')
    # #     # # ax1.set_xlabel('residence time [s]')
    # #     # # ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
    # #     # # ax2.set_ylabel('temperature [K]', color='C1')
    # #     f.tight_layout()
    # #     plt.show()






if __name__ == "__main__":
    run_combustor_atm_sim()
