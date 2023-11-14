import os
import sys
import click
import cantera as ct
import matplotlib.pyplot as plt
from cac.constants import DATA_DIR
from cac.reactors import AerosolSolution, AerosolReactor, AerosolConstPressureReactor

@click.command()
@click.option('--equiv_ratio', default=1.0, help='Equivalence ratio for the fuel')
@click.option('--test', default=False, help='Simply just test integration without recording')
@click.option('--steadystate/--no-steadystate', default=False, help='Integrate directly to steady state')
@click.option('--precon/--no-precon', default=True, help='Add preconditioning to the simulation')
@click.option('--endtime', default=1.0, help='Run the simulation to the set end time')
@click.option('--nsteps', default=10, help='The number of steps to take between records')
def run_combustor_atm_sim(equiv_ratio, test, steadystate, precon, endtime, nsteps):
    """ a real jet-A of intended average composition,
        POSF 4658

    Args:
        folder (_type_): _description_
        test (_type_): _description_
        endtime (_type_): _description_
    """
    # append appropriate directories
    sys.path.append(DATA_DIR)
    sys.path.append(os.path.join(DATA_DIR, "jetfuel"))

    # creation of combustor portion of simulation
    fuel_model = os.path.join(DATA_DIR, f"scombustor.yaml")
    fuel = AerosolSolution(fuel_model, name="gas", transport_model=None)

    # Create a Reservoir for the inlet, set to a methane/air mixture at a specified
    # equivalence ratio
    X_fuel = "NC10H22:0.4267, IC8H18: 0.3302, C7H8: 0.2431, H2S:0.001"
    fuel.TP = 1200, 4 * ct.one_atm
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, "H2O: 0.04, O2:0.2095, N2:0.7808")
    fuel.equilibrate("HP")

    # create inlet fuel tank
    inlet = ct.Reservoir(fuel)

    # create combustor
    # fuel.equilibrate("HP")
    combustor = AerosolReactor(fuel)
    combustor.volume = 2.36e-2

    # create atmosphere model
    atms = AerosolSolution(fuel_model, name="gas", transport_model=None)
    atms.TPX = 300, ct.one_atm, "H2O: 0.04, O2:0.2095, N2:0.7808"
    atms.equilibrate('HP')

    # create outlet of atmosphere
    far_field = ct.Reservoir(atms)

    # create atmosphere reactor
    atmosphere = AerosolReactor(atms)
    atmosphere.volume = 1e6
    # Flow between reactors
    ctres = 0.005
    def mass_flow_combustor(t):
        return combustor.mass / ctres

    # connect inlet to combustor
    inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mass_flow_combustor)

    # combustor to atmosphere
    outlet_mfc = ct.MassFlowController(combustor, atmosphere, mdot=mass_flow_combustor)

    # atmosphere to far field
    # outlet_far = ct.MassFlowController(atmosphere, far_field, mdot=atmos_to_ff)
    # outlet_far = ct.PressureController(atmosphere, far_field, primary=outlet_mfc, K=0.1)
    # setup reactor network
    net = ct.ReactorNet([combustor, atmosphere])
    net.derivative_settings = {"skip-falloff": True,
                               "skip-third-bodies": True,
                               "skip-coverage-dependence": True,
                               "skip-electrochemistry": True,
                               "skip-flow-devices": True,
                               "skip-walls": True}
    if precon:
        net.preconditioner = ct.AdaptivePreconditioner()

    def pressure(r):
        return r.T * ct.gas_constant * r.mass / r.thermo.mean_molecular_weight / r.volume

    def print_reactor_stats():
        print(f"Integrated to {net.time}..")
        print(combustor.mass, atmosphere.mass)
        print(combustor.T, atmosphere.T)
        print(combustor.volume, atmosphere.volume)
        print(combustor.density, atmosphere.density)
        print(pressure(combustor), pressure(atmosphere))
        print()

    if steadystate:
        net.advance_to_steady_state()
        print(net.time)
    elif not test:
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        comb_states = ct.SolutionArray(fuel)
        atms_states = ct.SolutionArray(atms)
        times = []
        # adding the initial conditions
        comb_states.append(combustor.thermo.state)
        atms_states.append( atmosphere.thermo.state)
        times.append(net.time)
        # loop for an hour of simulation time
        while net.time < endtime:
            print_reactor_stats()
            # try:
            for i in range(nsteps):
                net.step()
            comb_states.append(combustor.thermo.state)
            atms_states.append(atmosphere.thermo.state)
            times.append(net.time)
            # except Exception as e:
            #     break
        print_reactor_stats()
    else:
        net.step()
    # investigations
    # C7H8, NA, SA
    # some plotting
    if not test:
        # Plot results
        f, ax1 = plt.subplots(1, 1)
        ax1.plot(times, comb_states("CH4").Y*combustor.mass, '.-', color='b')
        ax2 = ax1.twinx()
        ax2.plot(times, atms_states("O3").Y*atmosphere.mass, '.-', color='g')
        ax2.plot(times, atms_states("CO2").Y*atmosphere.mass, '.-', color='r')
        # ax2.plot(times, atms_states("NO2").Y, '--', color='g')
        # ax2.plot(times, atms_states("NO3").Y, '.-', color='g')
        # ax2.plot(times, atms_states.T, '.-', color='C1')
        # # ax1.set_xlabel('residence time [s]')
        # # ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
        # # ax2.set_ylabel('temperature [K]', color='C1')
        f.tight_layout()
        plt.show()


if __name__ == "__main__":
    run_combustor_atm_sim()
