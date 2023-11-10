import os
import sys
import click
import cantera as ct
import matplotlib.pyplot as plt
from cac.constants import DATA_DIR
from cac.reactors import AerosolSolution, AerosolReactor, AerosolConstPressureReactor

@click.command()
@click.argument('folder', nargs=1)
@click.option('--test', default=False, help='Number of greetings.')
@click.option('--steadystate/--no-steadystate', default=False, help='Number of greetings.')
@click.option('--endtime', default=1.0, help='Number of greetings.')
def run_combustor_atm_sim(folder, test, steadystate, endtime):
    """ a real jet-A of intended average composition,
        POSF 4658

    Args:
        folder (_type_): _description_
        test (_type_): _description_
        endtime (_type_): _description_
    """
    # append appropriate directories
    sys.path.append(os.path.abspath(folder))
    sys.path.append(DATA_DIR)
    # constants used by both systems
    entrainment_ratio = 0.1 # amount of fresh air flowing into atmosphere compared to combustor mass flow
    # creation of combustor portion of simulation
    fuel_model = os.path.join(DATA_DIR, f"combustor.yaml")
    fuel = AerosolSolution(fuel_model, name="jetfuel", transport_model=None)

    # Create a Reservoir for the inlet, set to a methane/air mixture at a specified
    # equivalence ratio
    equiv_ratio = 1.0  # stoichiometric combustion
    X_fuel = "NC10H22:0.4267, IC8H18: 0.3302, C7H8: 0.2431"
    fuel.TP = 1200, 2 * ct.one_atm
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, 'O2:1.0, N2:3.76')
    # create inlet fuel tank
    inlet = ct.Reservoir(fuel)
    outlet = ct.Reservoir(fuel)

    # create combustor
    # fuel.equilibrate("HP")
    combustor = AerosolReactor(fuel)
    combustor.volume = 1.0

    # create atmosphere model
    atms = AerosolSolution(fuel_model, name="jetfuel", transport_model=None)
    atms.TPX = 1200, ct.one_atm, "O2:1.0, N2:3.76"
    atms.equilibrate('HP')

    # create outlet of atmosphere
    far_field = ct.Reservoir(atms)
    entrainment = ct.Reservoir(atms)

    # create atmosphere reactor
    atmosphere = AerosolConstPressureReactor(atms)
    atmosphere.volume = 1.0
    # Flow between reactors
    ctres = 0.1
    atres = 0.5
    def mass_flow_combustor(t):
        return combustor.mass / ctres

    def mass_flow_atmosphere(t):
        return atmosphere.mass / atres

    # connect inlet to combustor
    inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mass_flow_combustor)

    # combustor to atmosphere
    outlet_mfc = ct.MassFlowController(combustor, atmosphere, mdot=mass_flow_combustor)
    # outlet_mfc = ct.PressureController(combustor, atmosphere, primary=inlet_mfc, K=0)

    # # entrainment to atmosphere
    # entrain_mfc = ct.MassFlowController(entrainment, atmosphere, mdot=entrainment_mdot)

    # # atmosphere to far field
    outlet_far = ct.PressureController(atmosphere, far_field, primary=outlet_mfc, K=0)
    # outlet_far = ct.PressureController(atmosphere, far_field, primary=entrain_mfc, K=0.01)

    # setup reactor network
    net = ct.ReactorNet([combustor, atmosphere])
    net.derivative_settings = {"skip-falloff": True,
                               "skip-third-bodies": True,
                               "skip-coverage-dependence": True,
                               "skip-electrochemistry": True,
                               "skip-flow-devices": True,
                               "skip-walls": True}
    net.preconditioner = ct.AdaptivePreconditioner()

    if steadystate:
        net.advance_to_steady_state()
        print(net.time)
    elif not test:
        # Run a loop over decreasing residence times, until the reactor is extinguished,
        # saving the state after each iteration.
        comb_states = ct.SolutionArray(fuel)
        atms_states = ct.SolutionArray(atms)
        times = []
        # loop for an hour of simulation time
        while net.time < endtime:
            print(f"Integrated to {net.time}..")
            print(combustor.mass, atmosphere.mass)
            comb_states.append(combustor.thermo.state)
            atms_states.append( atmosphere.thermo.state)
            times.append(net.time)
            net.step()
        # adding the last
        comb_states.append(combustor.thermo.state)
        atms_states.append( atmosphere.thermo.state)
        times.append(net.time)
    else:
        net.step()

    # some plotting
    if not test:
        # Plot results
        f, ax1 = plt.subplots(1, 1)
        ax1.plot(times, comb_states("CH4").Y*combustor.mass, '.-', color='b')
        ax2 = ax1.twinx()
        ax2.plot(times, atms_states("CH4").Y*atmosphere.mass, '.-', color='g')
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
