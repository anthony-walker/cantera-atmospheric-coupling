import os
import csv
import cantera as ct
import numpy as np
from cac.constants import DATA_DIR

class PlumeSolution(ct.Solution):
    "Wrapper to allow assignment of custom attributes"

class PlumeReactor(ct.ExtensibleIdealGasConstPressureMoleReactor):

    def __init__(self, label, *args, **kwargs):
        super(PlumeReactor, self).__init__(label, *args, **kwargs)
        self.state_air = []
        self.enthalpy = []
        self.cp = []
        self.er_start = 0 # starting index for omega since time only moves forward
        self.thermo.zenith_angle = 0
        # open entrainment rate data
        self.entrainment = True
        with open(os.path.join(DATA_DIR, "entrainment-rates.csv"), "r") as f:
            rdr = csv.reader(f, delimiter=",")
            next(rdr)
            next(rdr)
            data = zip(*[r for r in rdr])
            self.tX, self.wX, self.tT, self.wT = [list(map(float, filter(lambda x: x, l))) for l in data]

    @property
    def TX_air(self):
        return self.state_air

    @TX_air.setter
    def TX_air(self, value):
        self.state_air = np.array(value)


    @property
    def enthalpy_air(self):
        return self.enthalpy

    @enthalpy_air.setter
    def enthalpy_air(self, value):
        self.enthalpy = np.array(value)

    @property
    def cp_air(self):
        return self.cp

    @cp_air.setter
    def cp_air(self, value):
        self.cp = np.array(value)

    def log_rate_interpolation(self, x, x1, y1, x2, y2):
        nlog = np.log10
        logy = nlog(x/x1) * nlog(y2/y1) / nlog(x2/x1) + nlog(y1)
        return 10 ** logy

    def after_eval(self, t, LHS, RHS):
        # updating zenith angle after eval
        self.thermo.zenith_angle = np.clip(np.mod(2*np.pi*t / (24*60*60), 2*np.pi), 0, np.pi) - np.pi/2
        # entrainment model
        if self.entrainment:
            if len(self.state_air) == 0:
                raise Exception("PlumeReactor: TX_air not set")
            elif len(self.enthalpy) == 0:
                raise Exception("PlumeReactor: enthalpy_air not set")
            elif len(self.cp) == 0:
                raise Exception("PlumeReactor: cp_air not set")
            # calculate entrainment rate by logarithmic interpolation
            def entrain_loop(times, rates):
                for i in range(len(times) - 1):
                    if t < times[0]:
                        omega_ent = rates[0]
                        break
                    elif t >= times[i] and t < times[i+1]:
                        omega_ent = self.log_rate_interpolation(t, times[i], rates[i], times[i+1], rates[i+1])
                        break
                    elif t > times[-1]:
                        omega_ent = 1 / t
                return omega_ent
            omega_T = entrain_loop(self.tT, self.wT)
            omega_X = entrain_loop(self.tX, self.wX)
            # find most from ideal gas law
            T,P = self.thermo.TP
            Nt = P * self.volume / ct.gas_constant / self.T
            moles = self.state_air[1:] * Nt
            # addition of value to energy equation for thermal entrainment
            RHS[0] -= omega_T * (self.T - self.state_air[0]) * self.mass * self.thermo.cp_mass
            # addition of value for species equation entrainment
            for i in range(len(moles)):
                RHS[i+1] += omega_X * moles[i]


class DilutionReactor(ct.ExtensibleIdealGasConstPressureMoleReactor):
    def __init__(self, *args, mdot, beta_da, beta_mixing, **kwargs):
        # pop variables
        self.zloc = kwargs.pop("zloc", 0.0)
        self.mixing_scale = kwargs.pop("mixing_scale", 0.55)
        self.dilution_scale_start = kwargs.pop("dilution_scale_start", 0.95)
        self.dilution_scale_end = kwargs.pop("dilution_scale_end", 1.0)
        self.total_length = kwargs.pop("total_length", 0.075) # meters
        self.flow_area = kwargs.pop("flow_area", 0.15) # meters squared
        # initialize super class
        super().__init__(*args, **kwargs)
        # new variables
        self.mass_flow_rate = mdot
        # constants
        self.beta_da = beta_da
        self.beta_mixing = beta_mixing

    def get_beta(self, z):
        if self.zloc >= 0 and self.zloc <= self.total_length * self.mixing_scale:
            return self.beta_mixing
        elif self.zloc >= self.dilution_scale_start * self.total_length and self.zloc <= self.total_length * self.dilution_scale_end:
            return self.beta_da
        return 0

    def after_initialize(self, t0):
        self.n_vars += 2
        self.i_z = self.n_vars - 1
        self.i_mdot = self.n_vars - 2

    def after_get_state(self, y):
        y[self.i_z] = self.zloc
        y[self.i_mdot] = self.mass_flow_rate

    def after_update_state(self, y):
        self.zloc = y[self.i_z]
        self.mass_flow_rate = y[self.i_mdot]

    def after_eval(self, t, LHS, RHS):
        # dz/dt evaluation
        RHS[self.i_z] = self.mass_flow_rate / self.thermo.density / self.flow_area
        RHS[self.i_mdot] = self.get_beta(self.zloc)

    def before_component_index(self, name):
        if name == 'zloc':
            return self.i_z
        elif name == 'mdot':
            return self.mass_flow_rate

    def before_component_name(self, i):
        if i == self.i_z:
            return 'zloc'
        elif i == self.i_mdot:
            return 'mdot'

    def get_mass_flow_rate(self, t):
        return self.mass_flow_rate
