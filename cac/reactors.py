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
        if os.path.isfile("aero-ro2-sum-test.txt"):
            with open("aero-ro2-sum-test.txt") as f:
                content = f.read()
                self.ro2_species = [sp.strip() for sp in content.split("\n")[:-1]]
        else:
            self.ro2_species = []
        self.state_air = []
        self.times = []
        self.erates = []
        self.er_start = 0 # starting index for omega since time only moves forward
        # open entrainment rate data
        self.entrainment = True
        with open(os.path.join(DATA_DIR, "entrainment_rates.csv"), "r") as f:
            csv_data = csv.reader(f, delimiter=",")
            next(csv_data) # skip name row
            for t, r in csv_data:
                self.times.append(float(t))
                self.erates.append(float(r))

    @property
    def TX_air(self):
        return self.state_air

    @TX_air.setter
    def TX_air(self, value):
        self.state_air = np.array(value)

    def log_rate_interpolation(self, x, x1, y1, x2, y2):
        nlog = np.log10
        # logy = (nlog(x) - nlog(x1)) * (nlog(y2) - nlog(y1)) / (nlog(x2) - nlog(x1)) + nlog(y1)
        logy = nlog(x/x1) * nlog(y2/y1) / nlog(x2/x1) + nlog(y1)
        return 10 ** logy

    def before_eval(self, t, LHS, RHS):
        self.thermo.zenith_angle = np.clip(np.mod(2*np.pi*t / (24*60*60), 2*np.pi), 0, np.pi) - np.pi/2
        # determine the RO2 sum
        self.thermo.ro2_sum = 0
        for sp in self.ro2_species:
            self.thermo.ro2_sum += self.thermo.concentrations[self.thermo.species_index(sp)]
        if self.entrainment:
            if len(self.state_air) == 0:
                raise Exception("PlumeReactor: TX_air not set")
            # calculate entrainment rate by logarithmic interpolation
            for i in range(len(self.times) - 1):
                if t < self.times[0]:
                    omega_ent = self.erates[0]
                    break
                elif t >= self.times[i] and t < self.times[i+1]:
                    omega_ent = self.log_rate_interpolation(t, self.times[i], self.rates[i], self.times[i+1], self.rates[i+1])
                    break
                elif t > self.times[-1]:
                    omega_ent = 1 / t
            # find most from ideal gas law
            T,P = self.thermo.TP
            Nt = P * self.volume / ct.gas_constant / T
            moles = self.state_air * Nt
            print(omega_ent)
            # addition of value to energy equation for thermal entrainment
            RHS[0] -= omega_ent * (self.T - self.state_air[0])
            # addition of value for species equation entrainment
            for i in range(1, len(self.state_air)):
                RHS[i] += omega_ent * moles[i]
