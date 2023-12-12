import os
import re
import warnings
import requests
import string
import random
import ruamel.yaml
from rdkit import Chem
import pubchempy as pcp
import concurrent.futures as cf
from cac.constants import DATA_DIR
from cac.extractor.msearch import msearch
from rdkit.Chem import rdMolDescriptors
from cantera import ck2yaml

yaml = ruamel.yaml.YAML(typ='safe')
yaml.default_flow_style=False
MCM_SPECIES_URL = "https://mcm.york.ac.uk/MCM/species/{:s}"

# add species that are assumed available in atmospherics
add_species = ["H2O", "O2", "H", "N2", "AR", "CO2", "HCL", "HBR"]

# global variable species_data place holder
species_data = {}


with open(os.path.join(DATA_DIR, "functional-groups.yaml")) as f:
    functional_groups = yaml.load(f)["smiles-groups"]
functional_group_list = sorted(functional_groups.items(), key=lambda x: len(x[0]), reverse=True)

def generate_random_string(length=10):
    characters = string.ascii_letters + string.digits
    return ''.join(random.choice(characters) for _ in range(length))


def print_functional_group_formulas():
    with open(os.path.join(DATA_DIR, "functional-groups.yaml")) as f:
        fct_gs = yaml.load(f)["functional-groups"]
    for group, smiles in fct_gs.items():
        mol = Chem.MolFromSmiles(smiles)
        print(group, smiles, rdMolDescriptors.CalcMolFormula(mol))

def smiles_to_semi_structural(smiles):
    # Start by replacing functional groups
    groups = []
    temp_smiles = smiles.strip()
    for group_smiles, notation in functional_group_list:
        if group_smiles in temp_smiles:
            groups.append(notation)
            temp_smiles = temp_smiles.replace(group_smiles, "")
    if temp_smiles:
        return generate_random_string(length=20)
    condensed = "".join(groups)
    return condensed

def inchi_to_composition(inchi_string):
    mol = Chem.MolFromInchi(inchi_string)
    if not mol:
        raise ValueError("Invalid INCHI string.")
    # Initialize a dictionary to hold atom counts
    composition = {}
    # Iterate over atoms and count each type
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        composition[atom_symbol] = composition.get(atom_symbol, 0) + 1
        if atom.GetTotalNumHs() > 0:
            composition["H"] = composition.get("H", 0) + atom.GetTotalNumHs()
    return composition


def smiles_to_composition(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        raise ValueError("Invalid SMILES string.")
    # Initialize a dictionary to hold atom counts
    composition = {}
    # Iterate over atoms and count each type
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        composition[atom_symbol] = composition.get(atom_symbol, 0) + 1
        if atom.GetTotalNumHs() > 0:
            composition["H"] = composition.get("H", 0) + atom.GetTotalNumHs()
    return composition


def make_species_database(dir_name):
    files = os.listdir(dir_name)
    all_species = {}
    for cf in files:
        with open(os.path.join(dir_name, cf), "r") as f:
            data = yaml.load(f)
            if "species" in data.keys():
                for sp in data["species"]:
                    # ensure all names are upper cases
                    # yaml reads species NO as false
                    sp_name = sp["name"].upper() if sp["name"] else "NO"
                    sp["name"] = sp_name
                    # check if the species exists
                    if sp_name not in all_species:
                        all_species[sp_name] = sp
                    else:
                        all_species[sp_name].update(sp)
    # output species to file
    with open("all-model-species.yaml", "w") as f:
        yaml.dump(all_species, f)


def get_not_found_species(specie):
    # request html data from MCM
    # print(f"Trying to find {specie}...")
    found = False
    sp_data = species_data.get(specie, {})
    if sp_data:
        print(f"Found {specie}...")
        return sp_data
    else:
        construct_data = {"name": specie, "composition":{}, "thermo":{}}
        note_str = f"""{specie} was not found in all-species-data and was constructed as well as possible using web requests and assuming ideal gas thermo.\n"""
        # assume ideal gas thermo
        construct_data["thermo"] = {"model":"constant-cp"}
        response = requests.get(MCM_SPECIES_URL.format(specie))
        if response.status_code == 200:
            inchi = ""
            smiles = ""
            molecule = None
            # try to locate from smiles
            # get smiles string
            res = re.search(r'strong[\>]Smiles: [^\<]*', response.text)
            smiles = res.group(0)[15:].strip()
            inchi = res.group(0)[14:].strip()
            # try to get from rmg
            try:
                chemkin = msearch(smiles)
                assert chemkin != ""
                sp_data = convert_chemkin(specie, chemkin)
                print(f"RMG Found {sp_data['alias']}:{specie}")
                sp_data["smiles"] = smiles
                sp_data["inchi"] = inchi
                return sp_data
            except:
                pass
            # try to get smiles from group
            try:
                cmp = pcp.get_compounds(smiles, "smiles")[0]
                molecule = Chem.MolFromSmiles(smiles)
                ffstr = "SMILES"
            except:
                pass
            # try to get smiles from inchi
            res = re.search(r'strong[\>]InChI: InChI=[^\<]*', response.text)
            try:
                cmp = pcp.get_compounds(inchi, "inchi")[0]
                if not smiles:
                    smiles = cmp.canonical_smiles
                if molecule is None:
                    molecule = Chem.MolFromInchi(inchi)
                    ffstr = "INCHI"
            except:
                pass
            # get molecule name
            if molecule is not None:
                molecule_name = rdMolDescriptors.CalcMolFormula(molecule)
            else:
                molecule_name = generate_random_string(length=20)
            # get semi-structural name
            if smiles:
                semi_struct_name = smiles_to_semi_structural(smiles)
            else:
                semi_struct_name = generate_random_string(length=20)
            # check for the molecule name otherwise try the semi-structural
            if species_data.get(molecule_name, {}):
                sm_data = species_data[molecule_name]
                test_comp = smiles_to_composition(smiles)
            elif species_data.get(semi_struct_name, {}):
                ffstr = "SEMI-STRUCTURAL"
                molecule_name = semi_struct_name
                sm_data = species_data[semi_struct_name]
                test_comp = smiles_to_composition(smiles)
            else:
                sm_data = {}
                test_comp = {}
            # check the composition
            if sm_data:
                if test_comp == sm_data.get("composition", {}):
                    print(f"Found with {ffstr} {specie}...")
                    construct_data.update(sm_data)
                    construct_data["name"] = specie
                    construct_data["alias"] = molecule_name
                    note_str = ((construct_data.get("note", " ") + note_str).strip()
                                    + "\n")
                    note_str += f"It was found using the {ffstr} representation: {molecule_name}."
                    found = True
            # add smiles and inchi string to notes
            note_str += f"SMILES: {smiles}\n"
            note_str += f"InCHI: {inchi}\n"
            # try to construct composition if it doesn't have one
            if not construct_data.get("composition", {}):
                comp = {}
                try:
                    comp = inchi_to_composition(inchi)
                except:
                    pass
                if not comp:
                    try:
                        comp = smiles_to_composition(smiles)
                    except:
                        pass
                if not comp:
                    note_str += f"Composition not constructed as an exception was encountered: {specie}."
                    formula = input(f"Composition not found for {specie}, please enter formula:")
                # set composition
                construct_data["composition"] = comp
            construct_data["note"] = note_str
            if not found:
                print(f"Constructed as well as possible {specie}...")
            if smiles:
                construct_data["smiles"] = smiles
            if inchi:
                construct_data["inchi"] = inchi
        return construct_data


def assign_global_species_data():
    global species_data
    all_sp_file = os.path.join(DATA_DIR, "all-model-species.yaml")
    with open(all_sp_file, "r") as f:
        species_data = yaml.load(f)


def convert_chemkin(species, chemkin, cleanup=True, add_alias=True):
    # write
    ckf = f"{species}.txt"
    with open(ckf, "w") as f:
        f.write("THERMO ALL\n0 200.00  1000.00  6000.00\n")
        f.write(chemkin)
    # convert to yaml
    os.system(f"ck2yaml --thermo={ckf} --permissive")
    # read yaml
    yf = f"{species}.yaml"
    with open(yf, "r") as f:
        data = yaml.load(f)["species"][0]
    # remove files
    if cleanup:
        os.remove(ckf)
        os.remove(yf)
    # add alias
    if species and add_alias:
        data["alias"] = data["name"]
        data["name"] = species
    # return data
    return data


def get_species_data(prefix, dirname):
    sfile = os.path.join(dirname, f"{prefix}-species.txt")
    with open(sfile, "r") as f:
        species = f.read().split("\n")[:-1]
    species += add_species
    print("Loading all species data...")
    assign_global_species_data()

    print("Executing threads...")
    with cf.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        res = [executor.submit(get_not_found_species, s.upper()) for s in species]
        result_data = [r.result() for r in res]
    # store by name of species
    data = {}
    if len(species) != (len(result_data)):
        warnings.warn(f"""Species determined not the same size as species list, {len(species)} != {len(result_data)}""", UserWarning)
    for sp in result_data:
        data[sp["name"]] = sp
    return data

def write_species_extraction(prefix, dirname):
    data = get_species_data(prefix, dirname)
    # output species to file
    sfile = os.path.join(dirname, f"{prefix}-species.yaml")
    with open(sfile, "w") as f:
        yaml.dump(data, f)

if __name__ == "__main__":
    get_not_found_species("BZEMUCCO")
