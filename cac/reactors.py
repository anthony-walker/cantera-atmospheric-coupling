import os
import csv
import ephem
import datetime
import cantera as ct
import numpy as np
from cac.constants import DATA_DIR


class PlumeSolution(ct.Solution):
    "Wrapper to allow assignment of custom attributes"

class PlumeReactor(ct.ExtensibleIdealGasConstPressureMoleReactor):

    def __init__(self, label, *args, **kwargs):
        # get some possible kwargs
        self.latitude = kwargs.pop("latitude", 0)
        self.longitude = kwargs.pop("longitude", 0)
        self.start_day = kwargs.pop("start_day", 0)
        self.start_time = kwargs.pop("start_time", 0) # hours
        self.altitude = kwargs.pop("altitude", 0) # hours
        self.no_change_species = kwargs.pop("no_change_species", [])
        super(PlumeReactor, self).__init__(label, *args, **kwargs)
        self.state_air = []
        self.enthalpy = []
        self.cp = []
        self.er_start = 0 # starting index for omega since time only moves forward
        self.thermo.zenith_angle = self.calculate_solar_zenith_angle(0)
        self.zat = -1
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

    def calculate_solar_zenith_angle(self, t):
        # calculate day number from start day, start time, and time
        time = t / 3600 + self.start_time# convert to hours
        days_ahead = time // 24
        time = time % 24 # get time remainder
        day_number = self.start_day + days_ahead # calculate day number
        # Convert latitude and longitude to radians
        latitude_rad = np.radians(self.latitude)
        # Calculate the declination angle
        declination = 23.45 * np.sin(np.radians((360/365) * (day_number - 81)))
        # Calculate the time correction factor
        time_correction = 4 * (self.longitude - 15 * 1) + 60 * 1
        # Calculate the solar hour angle
        solar_hour_angle = 15 * (time - 12) + time_correction
        # Calculate the solar zenith angle
        solar_zenith_angle = np.arccos(np.sin(np.radians(latitude_rad)) * np.sin(np.radians(declination)) + np.cos(np.radians(latitude_rad)) * np.cos(np.radians(declination)) * np.cos(np.radians(solar_hour_angle)))
        # adjust sza for altitude
        altitude_angle = self.altitude / (6.371e6) # divide by earths radius
        solar_zenith_angle = np.arccos(np.sin(latitude_rad) * np.sin(solar_zenith_angle) + np.cos(latitude_rad) * np.cos(solar_zenith_angle) * np.cos(altitude_angle))
        solar_zenith_angle = np.clip(np.mod(solar_zenith_angle, 2*np.pi), 0, np.pi / 2)
        return solar_zenith_angle

    def pressure(self):
        return self.T * ct.gas_constant * self.mass / self.thermo.mean_molecular_weight / self.volume

    def replace_eval(self, t, LHS, RHS):
        # updating zenith angle after eval
        if t != self.zat:
            self.zat = t
            self.thermo.zenith_angle = self.calculate_solar_zenith_angle(t)
        # update M value
        self.thermo.M = self.pressure() / (self.T * ct.boltzmann) / (100**3) # molec / cm^3
        # evaluate original eval
        self.default_eval(t, LHS, RHS)
        # reduce change to 0 in no-change species
        if self.no_change_species:
            for sp in self.no_change_species:
                sid = self.component_index(sp)
                RHS[sid] = 0
                LHS[sid] = 1
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
        self.enthalpy = []
        self.mf_air = []
        self.enthalpy_mass = 0

    @property
    def partial_enthalpy_air(self):
        return self.enthalpy

    @partial_enthalpy_air.setter
    def partial_enthalpy_air(self, value):
        self.enthalpy = np.array(value)

    @property
    def enthalpy_mass_air(self):
        return self.enthalpy_mass

    @enthalpy_mass_air.setter
    def enthalpy_mass_air(self, value):
        self.enthalpy_mass = value

    @property
    def mass_fractions_air(self):
        return self.mf_air

    @mass_fractions_air.setter
    def mass_fractions_air(self, value):
        self.mf_air = np.array(value)

    def get_beta(self, z):
        if self.zloc >= 0 and self.zloc <= self.total_length * self.mixing_scale:
            return self.beta_mixing
        elif self.zloc >= self.dilution_scale_start * self.total_length and self.zloc <= self.total_length * self.dilution_scale_end:
            return self.beta_da
        return 0

    def after_initialize(self, t0):
        self.n_vars += 2
        self.i_z = self.n_vars - 2
        self.i_mdot = self.n_vars - 1

    def after_get_state(self, y):
        y[self.i_z] = self.zloc
        y[self.i_mdot] = self.mass_flow_rate

    def after_update_state(self, y):
        self.zloc = y[self.i_z]
        self.mass_flow_rate = y[self.i_mdot]

    def replace_eval(self, t, LHS, RHS):
        if len(self.partial_enthalpy_air) == 0:
            raise Exception("DilutionReactor: partial_enthalpy_air not set")
        elif self.enthalpy_mass_air == 0:
            raise Exception("DilutionReactor: enthalpy_mass_air not set")
        elif len(self.mass_fractions_air) == 0:
            raise Exception("DilutionReactor: mass_fractions_air not set")
        # current beta
        beta = self.get_beta(self.zloc)
        mws = self.thermo.molecular_weights
        Y = self.thermo.Y
        n = self.get_state()
        rates = self.kinetics.net_production_rates
        brhoa = beta / self.thermo.density / self.flow_area
        # added equations for dilution
        RHS[self.i_z] = self.mass_flow_rate / self.thermo.density / self.flow_area if self.zloc < self.total_length else 0
        RHS[self.i_mdot] = beta * RHS[self.i_z]
        # update to energy equation for dilution
        LHS[0] = self.thermo.cp_mass
        RHS[0] = - np.sum(self.thermo.partial_molar_enthalpies * rates) / self.thermo.density
        RHS[0] += brhoa * (self.enthalpy_mass_air - np.sum(self.thermo.partial_molar_enthalpies / self.thermo.molecular_weights * self.mass_fractions_air))
        # update species equations based for dilution
        for i in range(self.thermo.n_species):
            RHS[i+1] = (self.mass_flow_rate * rates[i] / self.mass) # update existing RHS term
            # add beta terms
            RHS[i+1] += n[i+1] * brhoa + brhoa * self.mass_flow_rate / mws[i] * (self.mass_fractions_air[i] - Y[i])

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
