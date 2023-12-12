import cantera as ct
from cac.extractor.msearch import msearch
from cac.extractor.species_extractor import convert_chemkin

def mole_fraction_calculator(far, ap, bp, lm):

    # create initial fuel setting
    X_fuel = {"N-C10H22":0.4267, "I-C8H18":0.3302, "C7H8":0.2431}
    # add H2S
    # hydrogen sulfide
    Xs = (10 / 1e6)
    rp = 1 - Xs - far - ap - bp - lm
    # adjust fuels for farnesane blend
    X_fuel = {k:v * rp for k,v in X_fuel.items()}
    X_fuel.update({"iC15H32": far, "H2S":Xs, "C10H16a": ap, "C10H16b":bp, "C10H16":lm})
    Xf = ", ".join([f"{n}:{v}" for n, v in X_fuel.items()])
    # make cantera solution
    spec = ct.Species.list_from_file("species.yaml")
    spec_gas = ct.Solution(thermo='ideal-gas', kinetics="gas", species=spec)
    spec_gas.set_equivalence_ratio(1.0, Xf, 'O2:1.0, N2:3.76')
    xfi = ""
    for i, mf in enumerate(spec_gas.TPX[2]):
        xfi += f"\"{spec_gas.species_names[i]}\":{mf:.3e},\n"
    print(xfi[:-2])

if __name__ == "__main__":
    mole_fraction_calculator(0.05, 0.05, 0.05, 0.05)
