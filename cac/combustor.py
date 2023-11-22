import os
import sys
import click
import numpy
import warnings
import cantera as ct
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from cac.constants import DATA_DIR
from cac.reactors import AerosolSolution, AerosolReactor, AerosolConstPressureReactor

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

@click.command()
@click.option('--equiv_ratio', default=1.0, help='Equivalence ratio for the fuel')
@click.option('--test', default=False, help='Simply just test integration without recording')
@click.option('--steadystate/--no-steadystate', default=False, help='Integrate directly to steady state')
@click.option('--precon/--no-precon', default=True, help='Add preconditioning to the simulation')
@click.option('--endtime', default=1.0, help='Run the simulation to the set end time')
@click.option('--nsteps', default=10, help='The number of steps to take between records')
@click.option("--farnesene", default=0.05, help="Farnesene blend percentage")
@click.option("--sulfur", default=0.001, help="Sulfur amount in mole fraction")
def run_combustor_atm_sim(equiv_ratio, test, steadystate, precon, endtime, nsteps, farnesene, sulfur):
    """ a real jet-A of intended average composition,
        POSF 4658

    Args:
        folder (_type_): _description_
        test (_type_): _description_
        endtime (_type_): _description_
    """
    # append appropriate directories
    sys.path.append(DATA_DIR)
    # creation of combustor portion of simulation
    fuel_model = os.path.join(DATA_DIR, f"combustor.yaml")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, name="combustor", transport_model=None)
    # create initial fuel setting
    X_fuel = {"N-C10H22":0.4267, "I-C8H18":0.3302, "C7H8":0.2431}
    # add H2S
    if sulfur > 0:
        X_fuel.update({"H2S":sulfur})
    # adjust fuels for farnesene blend
    X_fuel = {k:v * (1-farnesene) for k,v in X_fuel.items()}
    # add farnesene blend
    Xfarne = 1 - numpy.sum([v for k,v in X_fuel.items()])
    if Xfarne > 0:
        X_fuel.update({"iC15H32": Xfarne})
    X_fuel = ", ".join([f"{k}:{v}" for k, v in X_fuel.items()])
    X_air = "H2O: 0.04, O2:0.2095, N2:0.7808"
    # combustor geometrical parameters
    radius = 0.15
    height = 0.28
    # residence time in combustor
    ctres = 0.0063
    # design parameters
    EPR = 25 # engine pressure ratio
    EGT = 800 # K, exhaust gas temperature
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    fuel.TPX = T_atm, p_atm, X_air
    # isentropic compression
    p1 = p_atm * EPR
    T1 = T_atm * (p1 / p_atm) ** (ct.gas_constant / fuel.cp_mole)
    # set input fuel conditions
    fuel.TP = T1, p1
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, X_air)
    combustor = ct.IdealGasConstPressureMoleReactor(fuel)
    combustor.volume = numpy.pi * radius * radius * height
    # create a wall for heat transfer
    spark = ct.Wall(ct.Reservoir(fuel), combustor, Q=6.975e6)
    # setup reactor network
    net = ct.ReactorNet([combustor,])
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
    # advance to residence time
    net.advance(ctres)
    # assign state post combustor
    T2, p2 = combustor.thermo.TP
    print(T2, p2)
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
    print(M3, T3, p3)
    T4 = T3 * (1 + gmone / 2 * (M3 ** 2)) / (1 + gmone / 2 * (M4 ** 2))
    print(T4, p4)
    fuel.TPY = T4, p4, combustor.Y

    # # creating reservoir and atmosphere from
    # # create atmosphere model
    atms = AerosolSolution(fuel_model, name="gas", transport_model=None)
    atms.TPX = 300, ct.one_atm, X_air
    far_field = ct.Reservoir(atms)
    def entrainment(t):
        pass
    # create atmosphere reactor
    atms.TPY = T4, p4, combustor.Y
    atmosphere = AerosolConstPressureReactor(atms)
    atmosphere.volume = combustor.mass / fuel.density
    # create inlet reservoir for atmosphere
    exhaust = ct.Reservoir(atms)
    # setup mass flow from reservoir to atmosphere
    gc = 1 # 1 in SI system
    mdot = p4 * A4 * M4 * numpy.sqrt(gamma * gc / ct.gas_constant * T4)
    exhaust_mfc = ct.MassFlowController(exhaust, atmosphere, mdot=mdot)
    entrainment_mfc = ct.MassFlowController(exhaust, atmosphere, mdot=mdot)
    # setup reactor network
    net = ct.ReactorNet([atmosphere])
    net.derivative_settings = {"skip-falloff": True,
                            "skip-third-bodies": True,
                            "skip-coverage-dependence": True,
                            "skip-electrochemistry": True,
                            "skip-flow-devices": True,
                            "skip-walls": True}
    # setup preconditioner
    if precon:
        net.preconditioner = ct.AdaptivePreconditioner()
    net.step()
    # # Run a loop over decreasing residence times, until the reactor is extinguished,
    # # saving the state after each iteration.
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
