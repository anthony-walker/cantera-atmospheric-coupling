import os
import csv
import h5py
import numpy
import warnings
import cantera as ct
import multiprocessing as mp
import matplotlib.pyplot as plt
from cac.constants import DATA_DIR, COLORS
from cac.combustor import multizone_combustor
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


def combustor_design_params():
    # design parameters
    equiv_ratio = 1.77
    EPR = 25 # engine pressure ratio
    p_atm = 0.2 * ct.one_atm
    T_atm = 240 # K
    X_air = "O2:1.0, N2:3.76"
    X_fuel = ", ".join([f'{k}:{v}'for k, v in {"NC10H22":0.4267, "IC8H18":0.3302, "C7H8":0.2431}.items()])
    # create path to fuel model
    fuel_model = os.path.join(DATA_DIR, "jetfuel.yaml")
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
    try:
        print(f"Running thrust-level: {x:.3f}")
        combustor, ifm = multizone_combustor(fuel, x, equiv_ratio, X_fuel, X_air, n_pz=21, name_id=f"-verification-{x:0.1f}")
    except:
        print(f"Failed for Thrust-level: {x:.3f}")
        return None
    return (combustor.thermo.state, combustor.mass, ifm, x)


def combustor_verification():
    # run in parallel
    fuel, equiv_ratio, X_fuel, X_air = combustor_design_params()
    thrust_levels = numpy.linspace(0.1, 1.0, 10)
    # with mp.Pool(os.cpu_count()) as pool:
    #     res = pool.map(parallel_run_combustor, thrust_levels)
    # # solution array
    # combustor_states = ct.SolutionArray(fuel, extra=["mass", "TL", "phi", "mdot", "fuel_mass"])
    # for r in res:
    #     if r is not None:
    #         state, mass, ifm, tl = r
    #         combustor_states.append(state, mass=mass, TL=tl, phi=equiv_ratio, mdot=1.086*tl, fuel_mass=ifm)
    # # write out state and reformats
    cbs = f"combustor-verification.hdf5"
    # combustor_states.save(cbs, overwrite=True, name="thermo")
    # create plot with generated data
    with h5py.File(cbs, "r") as hf:
        Y_data = hf["thermo"]["data"]["Y"]
        mass_data = hf["thermo"]["data"]["mass"]
        fuel_mass_data = hf["thermo"]["data"]["fuel_mass"]
        mdot_data = hf["thermo"]["data"]["mdot"]
        thrust_levels = hf["thermo"]["data"]["TL"]
        # compute NOx
        mass_nox = Y_data[:, fuel.species_index("NO2")] * mass_data * 1000 # grams
        nox_ratio = mass_nox / fuel_mass_data
        print(numpy.array(fuel_mass_data))
        print(nox_ratio)

        print(mass_nox)
if __name__ == "__main__":
    mixing_box_model_verification()
