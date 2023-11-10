import re
import os
import click
import ruamel.yaml
import pubchempy as pcp
from rdkit import Chem
from cac.constants import DATA_DIR
from ruamel.yaml.comments import CommentedMap

def semi_structural_to_smiles(semi_structural, functional_group_list):
    # Start by replacing functional groups
    groups = []
    temp_struct = semi_structural.strip()
    for group_smiles, notation in functional_group_list:
        if notation in temp_struct:
            groups.append(group_smiles)
            temp_struct = temp_struct.replace(notation, group_smiles)
    smiles = "".join(groups)
    if temp_struct:
        raise Exception("Semi-structural formula not completed.")
    return smiles

def update_mechanism_with_smiles_and_inchi(filename):
    # open all data
    yaml = ruamel.yaml.YAML(typ='safe', pure=True)
    yaml.preserve_quotes = True
    yaml.sort_keys = False
    # get functional group data:
    with open(os.path.join(DATA_DIR, "functional-groups.yaml"), "r") as f:
        functional_groups = yaml.load(f)["smiles-groups"]
    functional_group_list = sorted(functional_groups.items(), key=lambda x: len(x[1]), reverse=True)
    with open(filename, "r") as f:
        mech_data = yaml.load(f)
    # get inchi and smiles data
    species = mech_data["species"]
    for sp in species:
        name = sp["name"]
        # check if they have smiles or inchi
        smiles = sp.get("smiles", "")
        inchi = sp.get("inchi", "")

        if name.startswith("I"):
            name = name.replace("I", "iso-")
        elif name.startswith("N"):
            name = name.replace("N", "n-")
        elif name.startswith("T"):
            name = name.replace("T", "t-")
        # try to get compound by name from pubchem
        print(f"******************{name}********************")

        if smiles:
            print(f"smiles exists {smiles}")
        if inchi:
            print(f"inchi exists {inchi}")
        try:
            if not smiles or not inchi:
                compound = pcp.get_compounds(name, "name")[0]
                if not smiles:
                    print(f"Found {smiles} from name via pubchem...")
                    smiles = compound.canonical_smiles
                if not inchi:
                    print(f"Found {inchi} from name via pubchem...")
                    inchi = compound.inchi
        except Exception as e:
            pass
        # # try to get smiles from semistructural
        # try:
        #     if not smiles:
        #         rsm = semi_structural_to_smiles(name, functional_group_list)
        #         mol = Chem.MolFromSmiles(rsm)
        #         smiles = rsm
        #         print(f"Found {smiles} from semi-structural name via rdkit...")
        #     if not inchi:
        #         inchi = Chem.MolToInchi(mol)
        #         print(f"Found {inchi} from semi-structural name via rdkit...")
        # except Exception as e:
        #     pass
        # try to get smiles and inchi from name next
        try:
            if not smiles or not inchi:
                mol = Chem.MolFromMolBlock(name)
                if not smiles:
                    smiles = Chem.MolToSmiles(mol)
                    print(f"Found {smiles} from name via rdkit...")
                if not inchi:
                    inchi = Chem.MolToInchi(mol)
                    print(f"Found {inchi} from name via rdkit...")
        except Exception as e:
            pass
        # # try to get smiles and inchi from composition if all else fails
        # try:
        #     if smiles and not inchi:
        #         compound = pcp.get_compounds(smiles, "smiles")[0]
        #     elif inchi and not smiles:
        #         compound = pcp.get_compounds(inchi, "inchi")[0]
        #     elif not smiles and not inchi:
        #         print("Found with composition")
        #         comp = sp["composition"]
        #         formula = ""
        #         for ele, ct in comp.items():
        #             formula += f"{ele}{ct}" if ct > 1 else ele
        #         compound = pcp.get_compounds(formula, "formula")[0]
        #     else:
        #         compound = None
        #     # set smiles and inchi
        #     if not smiles:
        #         smiles = compound.canonical_smiles
        #     if not inchi:
        #         inchi = compound.inchi
        # except Exception as e:
        #     pass
        # # add surface specification
        # if "(S)" in name:
        #     smiles = f"{smiles}-surface"
        #     inchi = f"{inchi}-surface"
        if not smiles:
            print(f"Smiles not found for {name}, input a value:")
        if not inchi:
            print(f"INCHI not found for {name}, input a value:")
        # add back to species
        if smiles:
            sp["smiles"] = smiles
        if inchi:
            sp["inchi"] = inchi
    updated = os.path.join(os.path.dirname(filename), f"updated-{os.path.basename(filename)}")
    with open(updated, "w") as f:
        yaml.dump(mech_data, f)


def spmatch(spone, sptwo):
    res1 = spone.get("smiles", "ONE") == sptwo.get("smiles", "TWO")
    res2 = spone.get("inchi", "ONE") == sptwo.get("inchi", "TWO")
    res3 = spone.get("name") == sptwo.get("name")
    return res1 or res2 or res3


def merge_mechanisms(fuel, atmosphere):
    """ A code to merge atmospheric and gas phase mechanisms

    Args:
        fuel (str): File name of fuel mechanism
        atmosphere (str): File name of atmosphere mechanism
    """
    yaml = ruamel.yaml.YAML(typ='safe', pure=True)
    # open all data
    with open(fuel, "r") as f:
        fuel_data = yaml.load(f)
    with open(atmosphere, "r") as f:
        atms_data = yaml.load(f)
    # get first phase for each set and merge them
    fphase = fuel_data["phases"][0]
    aphase = atms_data["phases"][1]
    # get merged elements
    merged_elements = list(set(fphase.get("elements", []) + aphase.get("elements", [])))
    # go through species and identify matching species
    matches = []
    all_matched = []
    for sp1 in fuel_data["species"]:
        for sp2 in atms_data["species"]:
            if spmatch(sp1, sp2):
                if sp1["name"] in all_matched or sp2["name"] in all_matched:
                    found = ()
                    for m1, m2 in matches:
                        if m1 == sp1["name"] or m2 == sp2["name"]:
                            found = (m1, m2)
                            break
                    if found:
                        if sp1["name"] == sp2["name"]:
                            matches.remove(found)
                            matches.append((sp1["name"], sp2["name"]))
                        elif found[0] == found[1]:
                            pass
                        else:
                            keep = input(f"Keep {sp1['name']}, {sp2['name']} or {found[0]}, {found[1]}: ")
                            if not keep:
                                matches.remove(found)
                                matches.append((sp1["name"], sp2["name"]))
                    else:
                        matches.append((sp1["name"], sp2["name"]))
                else:
                    matches.append((sp1["name"], sp2["name"]))
                all_matched.append(sp1["name"])
                all_matched.append(sp2["name"])
    # now make substitutions in atmospheric file and update phases etc
    matches = {asp:fsp for fsp, asp in matches}
    all_species = fuel_data["species"]
    for asp in atms_data["species"]:
        if matches.get(asp["name"], None) is None:
            all_species.append(asp)
    # species name list
    all_species_names = [sp["name"] for sp in all_species]
    # get all reactions
    all_reactions = fuel_data["reactions"]
    for r in atms_data["atmosphere-reactions"] + atms_data["aerosol-reactions"]:
        for asp, fsp in matches.items():
            eqn = f" {r['equation']} "
            if re.search(f"[ ]{asp}[ ]", eqn):
                res = re.search(f"[ ]{asp}[ ]", eqn)
                r["equation"] = re.sub(f"[ ]{asp}[ ]", f" {fsp} ", eqn).strip()
        # don't add duplicates
        found = False
        if not r.get("duplicate", False):
            for cr in all_reactions:
                eqn1 = cr["equation"]
                eqn2 = r["equation"]
                eqn1 = sorted(list(set(re.sub("[^A-Za-z]+", " ", eqn1).split())))
                eqn2 = sorted(list(set(re.sub("[^A-Za-z]+", " ", eqn2).split())))
                if eqn1 == eqn2:
                    found = True
                    break
        if not found:
            all_reactions.append(r)
    # update fuel phase
    fphase["species"] =  all_species_names
    fphase["elements"] = merged_elements
    fuel_data["extensions"] = atms_data["extensions"]
    fuel_data["species"] = all_species
    fuel_data["reactions"] = all_reactions
    fuel_data["phases"][0] =  fphase
    # add all info to file
    comb_file = fuel.split(".")[0] + "-"+ atmosphere.split(".")[0] +".yaml"
    with open(comb_file, "w") as f:
        yaml.dump(fuel_data, f)

@click.command()
@click.argument('fuel', nargs=1, type=str)
@click.argument('atmosphere', nargs=1, type=str)
def merge_commandline(fuel, atmosphere):
    merge_mechanisms(fuel, atmosphere)
