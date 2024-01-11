import re
import os
import sys
import inspect
import numpy as np
import cantera as ct
import importlib.util as iu
from cac.constants import DATA_DIR

mcm_complex_rates = None
mcm_complex_funcs = None


class ScaledArrheniusData(ct.ExtensibleRateData):
    __slots__ = ("M", "T")

    def __init__(self):
        self.M = None
        self.T = None

    def update(self, gas):
        if self.M != gas.M or self.T != gas.T:
            self.M = gas.M
            self.T = gas.T
            return True
        else:
            return False



@ct.extension(name="scaled-arrhenius", data=ScaledArrheniusData)
class ScaledArrheniusRate(ct.ExtensibleRate):
    __slots__ = ("A", "b", "Ea", "scalar", "conversion")

    def set_parameters(self, node, units):
        system = ct.UnitSystem(node.get("units"))
        self.A = node["A"]
        self.b = node["b"]
        self.Ea = node["Ea"]
        self.scalar = node.get("scalar", 1)
        self.conversion = system.convert_rate_coeff_to(1, units)

    def get_parameters(self, node):
        node["A"] = self.A
        node["b"] = self.b
        node["Ea"] = self.Ea
        node["scalar"] = self.scalar

    def eval(self, data):
        return self.A * np.exp(-self.Ea / data.T) * self.conversion


class ZenithAngleData(ct.ExtensibleRateData):
    __slots__ = ("zenith_angle", "cza")

    def __init__(self):
        self.zenith_angle = None

    def update(self, gas):
        if self.zenith_angle != gas.zenith_angle:
            self.zenith_angle = gas.zenith_angle
            self.cza = np.cos(self.zenith_angle)
            return True
        else:
            return False


@ct.extension(name="zenith-angle-rate", data=ZenithAngleData)
class ZenithAngleRate(ct.ExtensibleRate):
    __slots__ = ("ell", "m", "n", "scalar", "conversion")

    def set_parameters(self, node, units):
        self.ell = node["l"]
        self.m = node["m"]
        self.n = node["n"]
        self.scalar = node.get("scalar", 1)
        system = ct.UnitSystem(node.get("units"))
        self.conversion = system.convert_rate_coeff_to(1, units)

    def get_parameters(self, node):
        node["l"] = self.ell
        node["m"] = self.m
        node["n"] = self.n
        node["scalar"] = self.scalar

    def eval(self, data):
        rate = self.ell * data.cza**self.m * np.exp(-self.n / data.cza) * self.scalar * self.conversion
        return rate


class FunctionData(ct.ExtensibleRateData):
    __slots__ = ("thermo")

    def __init__(self):
        self.thermo = None

    def update(self, gas):
        if self.thermo != gas:
            self.thermo = gas
            return True
        else:
            return False

@ct.extension(name="function-rate", data=FunctionData)
class FunctionRate(ct.ExtensibleRate):
    __slots__ = ("function", "pyfile", "ro2_species", "fargs", "conversion", "mspecies")

    def set_parameters(self, node, units):
        fcn_name = node["function"]
        # get ro2_species
        ro2_sumfile = node.get("ro2file", "")
        prefix = ro2_sumfile.split("-")[0]
        ro2_sumfile = os.path.join(DATA_DIR, prefix, ro2_sumfile)
        with open(ro2_sumfile, "r") as f:
            content = f.read()
            self.ro2_species = [sp.strip() for sp in content.split("\n")[:-1]]
        # this accesses a global module and assigns functions
        global mcm_complex_rates
        global mcm_complex_funcs
        if mcm_complex_rates is None:
            pyfile = node.get("pyfile")
            prefix = pyfile.split("_")[0]
            filename = os.path.join(DATA_DIR, prefix, pyfile)
            modname = pyfile.split(".")[0]
            spec = iu.spec_from_file_location(modname, filename)
            mcm_complex_rates = iu.module_from_spec(spec)
            spec.loader.exec_module(mcm_complex_rates)
            mcm_complex_funcs = dict(inspect.getmembers(mcm_complex_rates, inspect.isfunction))
        # assign function
        self.function = mcm_complex_funcs[fcn_name]
        self.fargs = inspect.getargspec(self.function).args
        system = ct.UnitSystem(node.get("units"))
        self.conversion = system.convert_rate_coeff_to(1, units)


    def get_parameters(self, node):
        node["function"] = self.function.__name__


    def eval(self, data):
        # set ro2 sum to 0
        ro2_frac = 0
        if "RO2" in self.fargs:
            for sp in self.ro2_species:
                ro2_frac += data.thermo.Y[data.thermo.species_index(sp)]
        # evaluate the function
        built_args = []
        for a in self.fargs:
            if "T" == a:
                built_args.append(data.thermo.T)
            elif "RO2" == a:
                built_args.append(ro2_frac * data.thermo.M)
            elif  "M" == a:
                built_args.append(data.thermo.M)
            else:
                built_args.append(data.thermo.Y[data.thermo.species_index(a)] * data.thermo.M)
        # convert based on argument length
        rate = self.function(*built_args) * self.conversion
        return rate
