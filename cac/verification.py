import os
import csv
import numpy
import warnings
import cantera as ct
import matplotlib.pyplot as plt
from cac.constants import DATA_DIR, COLORS
from cac.combustor import braggs_combustor
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
    with open(os.path.join(DATA_DIR, "box-model-verification.csv"), "r") as f:
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
    plt.savefig("box-model-verification.pdf")


def combustor_verification():
    # design parameters
    equiv_ratio = 0.25
    EPR = 25 # engine pressure ratio
    EGT = 800 # K, exhaust gas temperatures
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    X_air = "O2:1.0, N2:3.76"
    X_fuel = "H2:1.0"
    # create path to fuel model
    fuel_model = "h2o2.yaml"
    # creation of fuel thermo object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fuel = ct.Solution(fuel_model, transport_model=None, basis="mole")
    fuel.TPX = T_atm, p_atm, X_air
    # isentropic compression
    p1 = p_atm * EPR
    T1 = T_atm * (p1 / p_atm) ** (ct.gas_constant / fuel.cp_mole)
    # setup fuel for braggs combustor
    # set input fuel conditions
    fuel.TP = T1, p1
    fuel.set_equivalence_ratio(equiv_ratio, X_fuel, X_air, basis="mole")
    # run braggs combustor
    combustor, combustor_states, wsr, wsr_states = braggs_combustor(fuel)
    # plots to verify combustor behavior
    Y = numpy.array(combustor_states.X)
    ch4_id = combustor.thermo.species_index("H2")
    # ch4_id = combustor.thermo.species_index("")
    # plot data
    f, ax1 = plt.subplots(1, 1)
    ax1.loglog(combustor_states.t, Y[:, ch4_id], color=COLORS[3])
    ax2 = ax1.twinx()
    ax2.semilogx(combustor_states.t, combustor_states.T, color=COLORS[1])
    # x axis
    ax1.set_xlabel("Time [s]")
    # ax1 y axis setup
    ax1.set_ylabel("$Y_{H2}$", color=COLORS[3])
    ax1.tick_params(axis="y", colors=COLORS[3])
    # ax2 y axis setup
    ax2.set_ylabel("Temperature [K]", color=COLORS[1])
    ax2.tick_params(axis="y", colors=COLORS[1])
    plt.savefig("h2o2-combustor.pdf")
    # plt.plot(combustor_states.z, combustor_states.T)
    plt.show()


if __name__ == "__main__":
    mixing_box_model_verification()
