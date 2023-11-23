import os
import cantera as ct
import numpy as np

class PlumeSolution(ct.Solution):
    "Wrapper to allow assignment of custom attributes"

class PlumeReactor(ct.ExtensibleIdealGasConstPressureMoleReactor):

    def __init__(self, label, *args, **kwargs):
        super(PlumeReactor, self).__init__(label, *args, **kwargs)
        if os.path.isfile("aero-ro2-sum-test.txt"):
            with open("aero-ro2-sum-test.txt") as f:
                content = f.read()
                self.ro2_species = [sp.strip() for sp in content.split("\n")[:-1]]
        else:
            self.ro2_species = []
        self.T_amb = kwargs.get("T_amb", 300) # ambient temperature in kelvin
        self.t1 = kwargs.get("t1", 8) # seconds from karcher box model
        self.t2 = kwargs.get("t2", 66) # seconds from karcher box model
        self.Cv = kwargs.get("Cv", 3) # from karcher box model
        self.m = kwargs.get("m", 2) # from karcher box model
        self.n = kwargs.get("n", 50) # from karcher box model
        self.omega_X = 0
        self.omega_T = 0


    def before_eval(self, t, LHS, RHS):
        # Sample version of this that has the right periodicity. A full implementation would use the
        # latitude and longitude. Clipping accounts for the period where the sun is behind the earth.
        self.thermo.zenith_angle = np.clip(np.mod(2*np.pi*t / (24*60*60), 2*np.pi), 0, np.pi) - np.pi/2
        # determine the RO2 sum
        self.thermo.ro2_sum = 0
        for sp in self.ro2_species:
            self.thermo.ro2_sum += self.thermo.concentrations[self.thermo.species_index(sp)]
        # calculate omega
        if (t > 0):
            self.omega_X = np.zeros(self.thermo.n_species)
            self.omega_T = 0
            # TODO: ADD NEEDED DERIVATIVES HERE
        if (t > self.t1):
            self.omega_X = self.omega_X / (1 + self.Cv * self.omega_X * (t - self.t1))
            self.omega_T = self.omega_T / (1 + self.Cv * self.omega_T * (t - self.t1))
        if (t > self.t2):
            Cd = 1 / (1 + self.m * np.exp(-(t - self.t2) / (self.n * self.t2)))
            assert Cd <= 1
            self.omega_X = self.omega_X / (1 + Cd * self.omega_X * (t - self.t2))
            self.omega_T = self.omega_T / (1 + Cd * self.omega_T * (t - self.t2))

        # addition of value to energy equation for thermal entrainment
        RHS[0] -= self.omega_T * (self.T -  self.T_amb)
        # addition of value for species equation entrainment
        # for i, oi in enumerate(self.omega_X):
        #     RHS[i] -= oi * (x_i)
