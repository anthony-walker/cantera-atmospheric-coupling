## Cantera Atmospheric Coupling Documentation
In this repository is all the code I have used to test and develop my couplings of atmospheric chemistry to Cantera.
All saved or stored results for studies can be found in the `study` directory.

### Getting Started
The easiest way to repeat or use my analyses is to install this package with

```sh
pip install -e .
```
from the github repository.

### Verifications
All verifications can be run with `verify_*`. The outputs for these commands default to `./cac/data/verification`.

- Box model: `verify_box_model`
- Atmospheric model `verify_atmosphere`
- Combustor model `verify_combustor`

### Tools
A set of tools was developed in the process to perform these analyses, namely `extractor` and `merger` commands.
- `extractor`: This converts a model extracted from the Master Chemical Mechanisms as a `*.fac` file to a `*.yaml` for use in Cantera.
- `merger`: This merges any number of models based on the `smiles` and `inchi` fields in each species element.

### Simulation
The setup simulation in its latest state can be run with the `combustor` command.

```sh
Usage: combustor [OPTIONS]

Options:
  --equiv_ratio FLOAT  Equivalence ratio for the fuel
  --nsteps INTEGER     The number of steps to take between records
  --farnesane FLOAT    Farnesane blend percentage
  --sulfur FLOAT       Sulfur amount in mole fraction
  --outdir TEXT        Output directory for data
  --help               Show this message and exit.
```
