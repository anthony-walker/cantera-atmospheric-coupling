import re
import os
import sys
import copy
import click
import signal
import datetime
import warnings
import ruamel.yaml
from ruamel.yaml.comments import CommentedMap, CommentedSeq
import pubchempy as pcp
from rdkit import Chem
import cantera as ct
from cac.constants import DATA_DIR
from cac.reactors import PlumeSolution
from cac.rates import *
from ruamel.yaml.comments import CommentedMap

def FlowMap(*args, **kwargs):
    m = CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def FlowList(*args, **kwargs):
    lst = CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst


def spmatch(spone, sptwo):
    res1 = spone.get("smiles", "ONE") == sptwo.get("smiles", "TWO")
    res2 = spone.get("inchi", "ONE") == sptwo.get("inchi", "TWO")
    res3 = spone.get("name") == sptwo.get("name")
    res4 = spone.get("name") == sptwo.get("alias", "TWO")
    res5 = spone.get("alias", "ONE") == sptwo.get("name", "TWO")
    res = res1 or res2 or res3 or res4 or res5
    res = res and (spone["composition"] == sptwo["composition"])
    return res


def merge_all_mechanisms(*mechanisms):
    yaml = ruamel.yaml.YAML()
    # yaml.preserve_quotes = True
    yaml.default_flow_style=False
    # yaml.sort_keys = False
    data_sets = [yaml.load(open(m, "r")) for m in mechanisms]
    mnames = [m.split(".")[0].strip() for m in mechanisms]
    # get a complete list of elements
    elements = []
    for ds in data_sets:
        for ele in ds["phases"][0]["elements"]:
            elements.append(ele)
    elements = list(set(elements))
    # go through species and identify matching species
    primary_species = data_sets[0]["species"]
    smaps = []
    for i, ds in enumerate(data_sets[1:]):
        sp_map = {}
        csp = []
        for sp2 in ds["species"]:
            found = False
            for sp1 in primary_species:
                if spmatch(sp1, sp2) and sp1["name"] not in csp:
                    found = True
                    print(f"{sp1['name']} and {sp2['name']} match")
                    # transfer any missing keys
                    for k in sp2.keys():
                        if k not in sp1.keys():
                            sp1[k] = sp2[k]
                    # if names do not match add a mapping
                    if sp1["name"] != sp2["name"]:
                        sp_map[sp1["name"]] = sp2["name"]
                    break
            if not found:
                primary_species.append(sp2)
                csp.append(sp2["name"])
        smaps.append(sp_map)
    # get all names
    all_species_names = [sp["name"] for sp in primary_species]
    # get all reactions
    reaction_keys = []
    for i, ds in enumerate(data_sets):
        rkeys = list(filter(lambda x: "reaction" in x, ds.keys()))
        new_rkeys = [ck if ck != "reactions" else f"{mnames[i]}-reactions" for ck in rkeys]
        for i in range(len(rkeys)):
            ds[new_rkeys[i]] = ds.pop(rkeys[i])
        reaction_keys.append(new_rkeys)
    # make all reactions dictionary
    all_reactions = {rk:data_sets[0][rk] for rk in reaction_keys[0]}
    # make all the necessary substitutions
    for i, mp, rkeys in zip(range(1, 1 + len(smaps)), smaps, reaction_keys[1:]):
        for ps, nps in mp.items():
            for r in rkeys:
                for j, rct in enumerate(data_sets[i][r]):
                    eqn = f" {rct['equation']} "
                    if re.search(f"[ ]{nps}[ ]", eqn):
                        rct["equation"] = re.sub(f"[ ]{nps}[ ]", f" {ps} ", eqn).strip()
    # make a list of all reactions
    for i, rkeys in enumerate(reaction_keys[1:], start=1):
        for rk in rkeys:
            all_reactions[rk] = data_sets[i][rk]
    # capture all extensions
    all_extensions = []
    for ds in data_sets:
        for ext in ds.get("extensions", []):
            all_extensions.append(ext)
    # update main set and write model
    main_data_set = data_sets[0]
    main_data_set["extensions"] = all_extensions
    main_data_set["species"] =  primary_species
    main_phase = main_data_set["phases"][0]
    main_phase["elements"] = FlowList(elements)
    main_phase["species"] = FlowList(all_species_names)
    for rk, reacts in all_reactions.items():
        main_data_set[rk] = reacts
    main_data_set["phases"][0] = main_phase
    made_str = f"Made from: {', '.join(mechanisms)}."
    desc = main_data_set.get("description", made_str)
    desc += f" {made_str}" if desc != made_str else ""
    main_data_set["description"] = desc
    # add all info to file
    with open("merged.yaml", "w") as f:
        yaml.dump({"description":main_data_set.pop("description", "")}, f)
        yaml.dump({"generator":main_data_set.pop("generator", "")}, f)
        yaml.dump({"input-files":main_data_set.pop("input-files", [])}, f)
        yaml.dump({"cantera-version":main_data_set.pop("cantera-version", "3.0")}, f)
        dt = main_data_set.pop("date", datetime.datetime.now())
        if not isinstance(dt, str):
            dt = dt.strftime("%m/%d/%Y, %H:%M:%S")
        yaml.dump({"date":dt}, f)
        yaml.dump({"units":main_data_set.pop("units", {})}, f)
        yaml.dump({"extensions":main_data_set.pop("extensions", "")}, f)
        yaml.dump({"phases":main_data_set.pop("phases", "")}, f)
        yaml.dump({"species":main_data_set.pop("species", "")}, f)
        for yk in main_data_set:
            yaml.dump({yk: main_data_set[yk]}, f)


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
    try:
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
        # closure to periodically write out
        def write_out():
            updated = os.path.join(os.path.dirname(filename), f"updated-{os.path.basename(filename)}")
            with open(updated, "w") as f:
                yaml.dump(mech_data, f)
        # register signal handler
        def crash_mech_handler(sig, frame):
            write_out()
            sys.exit(0)

        signal.signal(signal.SIGINT, crash_mech_handler)
        # get inchi and smiles data
        species = mech_data["species"]
        all_smiles = []
        all_inchi = []
        lsp = len(species)
        for spi, sp in enumerate(species):
            name = sp["name"]
            updated = False
            # check if they have smiles or inchi
            smiles = ""
            inchi = ""
            if not sp.get("smiles", "") and not sp.get("inchi", ""):
                searches = []
                if re.match(r"(I-|i|iso|i-|iso-|IC)", name):
                    searches.append("iso")
                if re.match(r"(n|N-|n-|NC)", name):
                    searches.append("n")
                if re.match(r"(t|T-|t-|TC)", name):
                    searches.append("tert")
                if re.match(r"(sec-|s-|S-)", name):
                    searches.append("sec")
                if re.match(r"(cy-|cy)", name):
                    searches.append("cyclo")
                # elif name.startswith("T"):
                #     name = name.replace("T", "t-")
                # try to get compound by name from pubchem
                print()
                print(f"****************** {spi}/{lsp} {name} ********************")
                print("Searches: ", searches)
                print("Composition: ", sp["composition"])
                print("Note: ", sp.get("note", ""))
                keep = 1
                sname = input("Input search name: ").strip()
                if sname == "break":
                    break
                while keep:
                    try:
                        print(f"Searching for {sname}")
                        cmps = pcp.get_compounds(sname, "formula")
                        print(f"Found {len(cmps)} compounds")
                        for ci, cmp in enumerate(cmps):
                            fcomp = {}
                            for a in cmp.atoms:
                                fcomp[a.element] = fcomp.get(a.element, 0) + 1
                            if fcomp == sp["composition"] and cmp.isomeric_smiles not in all_smiles and cmp.inchi not in all_inchi:
                                print()
                                print(f"****************** {spi}/{lsp} {name} ********************")
                                print("Searches: ", searches)
                                print("Composition: ", sp["composition"])
                                print("Note: ", sp.get("note", ""))
                                print(f"--------------Compound {ci}-------------------")
                                print(cmp.molecular_formula)
                                print(cmp.isomeric_smiles)
                                print(cmp.inchi)
                                print(cmp.iupac_name)
                                print(cmp.synonyms)
                                keep = int(input("Enter 0 to keep the entry: ").strip())
                            if not keep:
                                smiles = cmp.isomeric_smiles
                                inchi = cmp.inchi
                                break
                            if keep == 2:
                                keep = 1
                                break
                            elif keep == 3:
                                inchi = input("Enter inchi: ").strip()
                                smiles = input("Enter smiles: ").strip()
                                keep = 0
                                break
                        if not cmps:
                            raise Exception(f"{sname} not found.")
                        elif keep:
                            sname = input("Input a NEW search name: ").strip()
                    except Exception as e:
                        print(e)
                        keep = 1
                        sname = input("Search failed try again: ").strip()
                        if sname == "skip":
                            keep = 0
                            smiles = ""
                            inchi = ""
                        elif sname == "manual":
                            inchi = input("Enter inchi: ").strip()
                            smiles = input("Enter smiles: ").strip()
                            keep = 0
                            break
            else:
                print(f"InChi and Smiles exists for {name}.")
                print()
            # add back to species
            if smiles:
                all_smiles.append(smiles)
                sp["smiles"] = smiles
                updated = True
            if inchi:
                all_inchi.append(inchi)
                sp["inchi"] = inchi
                updated = True
            # write out after every entry to avoid loss during errors
            if updated:
                write_out()
    except Exception as e:
        print(e)
        write_out()
        sys.exit(0)

def print_reaction_props(r):
    print(r.reaction_type)
    print(r.products)
    print(r.reactants)
    print(dir(r.rate))
    print(r.rate.type)
    print(r.rate.sub_type)
    print(r.reversible)
    tb = r.third_body
    print(tb)
    if tb is not None:
        print(tb.efficiencies)
    print()

def reactmatch(rone, rtwo, kin, verb=False):
    if verb:
        print("Testing: ")
        print(rone["equation"])
        print(rtwo["equation"])
        print()

    ro1 = ct.Reaction.from_dict(rone, kin)
    ro2 = ct.Reaction.from_dict(rtwo, kin)
    # if both equations are thirdbody check more simply
    tb1 = ro1.third_body
    tb2 = ro2.third_body
    match = False
    if tb1 is not None and tb2 is not None:
        p1 = ro1.products
        p2 = ro2.products
        r1 = ro1.reactants
        r2 = ro2.reactants
        match = match or (p1 == p2 and r1 == r2) or (p1 == r2 and r1 == p2)
    # check for tb in both reactions
    if not match:
        prodone = []
        reactone = []
        if tb1 is None:
            prodone.append(ro1.products)
            reactone.append(ro1.reactants)
        else:
            for k in tb1.efficiencies.keys():
                cr = {k:1.0}
                cr.update(ro1.reactants)
                cp = {k:1.0}
                cp.update(ro1.products)
                prodone.append(cp)
                reactone.append(cr)
        prodtwo = []
        reacttwo = []
        if tb2 is None:
            prodtwo.append(ro2.products)
            reacttwo.append(ro2.reactants)
        else:
            for k in tb2.efficiencies.keys():
                cr = {k:1.0}
                cr.update(ro2.reactants)
                cp = {k:1.0}
                cp.update(ro2.products)
                prodtwo.append(cp)
                reacttwo.append(cr)
        for p1, r1 in zip(prodone, reactone):
            for p2, r2 in zip(prodtwo, reacttwo):
                match = match or (p1 == p2 and r1 == r2) or (p1 == r2 and r1 == p2)
                if match:
                    break
    if match:
        print("Matches:")
        print(rone["equation"])
        print(rtwo["equation"])
        print()
        # input()
    return match, rtwo.get("duplicate", False)


def update_duplicates_combustor_mechanism():
    yaml = ruamel.yaml.YAML()
    # yaml.preserve_quotes = True
    yaml.default_flow_style=False
    yaml.representer.ignore_aliases = lambda *args: True
    # yaml.sort_keys = False
    combustor_data = yaml.load(open("merged.yaml", "r"))
    cphase = combustor_data["phases"][0]
    cphase["name"] = "combustor"
    cphase["reactions"] = ["farnesane-reactions", "sulfur-reactions"]
    aphase = copy.deepcopy(cphase)
    aphase["name"] = "atmosphere"
    aphase["reactions"] = ["atmosphere-reactions", "aerosol-reactions", "photolysis-reactions"]
    combustor_data["phases"].append(aphase)
    cks = [rk for rk in combustor_data["phases"][0]["reactions"]]
    aks = [rk for rk in combustor_data["phases"][1]["reactions"]]
    # create solution object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = ct.Species.list_from_file("merged.yaml")
        spec_gas = PlumeSolution(thermo='ideal-gas', kinetics="gas", species=spec)
    # update combustor reactions
    combustor_reactions = []
    for ck in cks:
        for cr in combustor_data.get(ck, []):
            duplicate = False
            known = False
            for r in combustor_reactions:
                duplicate, known = reactmatch(r, cr, spec_gas)
                if duplicate and known:
                    r["duplicate"] = True
                    cr["duplicate"] = True
                    break
                elif duplicate:
                    break
            if (not duplicate) or known:
                combustor_reactions.append(cr)
    combustor_data["combustor-reactions"] = combustor_reactions
    cphase["reactions"] = ["combustor-reactions"]
    # update atmospheric reactions
    # atms_reactions = []
    # for ak in aks:
    #     for cr in combustor_data.get(ak, []):
    #         duplicate = False
    #         known = False
    #         for r in atms_reactions:
    #             duplicate, known = reactmatch(r, cr, spec_gas)
    #             if duplicate and known:
    #                 r["duplicate"] = True
    #                 cr["duplicate"] = True
    #                 break
    #             elif duplicate:
    #                 break
    #         if (not duplicate) or known:
    #             atms_reactions.append(cr)
    atms_reactions = [ r for ak in aks for r in combustor_data.get(ak, [])]
    combustor_data["atmosphere-reactions"] = atms_reactions
    aphase["reactions"] = ["atmosphere-reactions"]
    # pop other reactions
    for rtype in ["farnesane-reactions", "sulfur-reactions", "aerosol-reactions", "photolysis-reactions"]:
        combustor_data.pop(rtype, "")
    # write to file
    with open("combustor.yaml", "w") as f:
        yaml.dump(combustor_data, f)

@click.command()
@click.argument('mechanisms', nargs=-1)
def merge_commandline(mechanisms):
    merge_all_mechanisms(*mechanisms)
    update_duplicates_combustor_mechanism()

def test_rmatch():
    # create solution object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = ct.Species.list_from_file("merged.yaml")
        spec_gas = PlumeSolution(thermo='ideal-gas', kinetics="gas", species=spec)
    r0 = {"equation": "H + O2 <=> O + OH",
          "rate-constant": {"A": 4.577e+19, "b": -1.4, "Ea": 1.044e+05}}
    r1 = {"equation": "C2H2 (+M) <=> C2H2 (+M)",
          "rate-constant": {"A": 4.577e+19, "b": -1.4, "Ea": 1.044e+05}}
    r2 = {"equation": "C2H4 + CH2a <=> C2H3 + CH3",
          "type": "pressure-dependent-Arrhenius",
          "rate-constants": [{"A":  1.77e+19, "a":  6787.0, "P":  "0.01 atm", "b":  -1.95}, {"A":  1.68e+19, "a":  4310.0, "P":  "0.1 atm", "b":  -1.8}, {"A":  4.16e+24, "a":  9759.0, "P":  "1.0 atm", "b":  -3.19}, {"A":  7.89e+24, "a":  13894.0, "P":  "10.0 atm", "b":  -3.08}, {"A":  7.36e+29, "a":  23849.0, "P":  "100.0 atm", "b":  -4.28}], "duplicate": True}
    r3 = {"equation": "C2H4 + CH2a <=> C2H3 + CH3",
          "type": "pressure-dependent-Arrhenius",
          "rate-constants": [{"A":  4300000000000.0, "a":  -110.0, "P":  "0.01 atm", "b":  0.19}, {"A":  226000000000.0, "a":  48.0, "P":  "0.1 atm", "b":  0.54}, {"A":  4920000000.0, "a":  600.0, "P":  "1.0 atm", "b":  1.02}, {"A":  147000000.0, "a":  1228.0, "P":  "10.0 atm", "b":  1.33}, {"A":  81100000000.0, "a":  5507.0, "P":  "100.0 atm", "b":  0.55}], "duplicate": True}
    r4 = {"equation": "C2H2 + OH <=> CHCHOH", "rate-constant": {"A": 4.577e+19, "b": -1.4, "Ea": 1.044e+05}}#"rate-constants":[{"A": 2.9e+64, "Ea": 10009.0, "P": "0.01 atm", "b": -18.57}, {"A": 2.6e+33, "Ea": 6392.0, "P": "0.01 atm", "b": -7.36}, {"A": 4.7e+59, "Ea": 9087.0, "P": "0.025 atm", "b": -16.87}, {"A": 4.4e+32, "Ea": 5933.0, "P": "0.025 atm", "b": -7.02}, {"A": 1.2e+28, "Ea": 3724.0, "P": "0.1 atm", "b": -5.56}, {"A": 6.4e+42, "Ea": 11737.0, "P": "0.1 atm", "b": -9.96}, {"A": 1.9e+44, "Ea": 6299.0, "P": "1.0 atm", "b": -11.38}, {"A": 3.5e+31, "Ea": 6635.0, "P": "1.0 atm", "b": -6.2}, {"A": 1.5e+24, "Ea": 3261.0, "P": "10.0 atm", "b": -4.06}, {"A": 4.5e+31, "Ea": 8761.0, "P": "10.0 atm", "b": -5.92}, {"A": 6.2e+20, "Ea": 2831.0, "P": "100.0 atm", "b": -2.8}, {"A": 1.6e+29, "Ea": 9734.0, "P": "100.0 atm", "b": -4.91}]}
    print(reactmatch(r0, r4, spec_gas))
    # print(reactmatch(r2, r3, spec_gas))

def check_all_smiles_and_inchi_unique(filename):
    yaml = ruamel.yaml.YAML()
    yaml.default_flow_style=False
    with open(filename, "r") as f:
        mech_data = yaml.load(f)
    species = mech_data["species"]
    for i, sp0 in enumerate(species):
        for j, sp1 in enumerate(species):
            if i != j:
                if sp0.get("smiles","ZERO") == sp1.get("smiles","ONE"):
                    print(f"{sp0['name']}: {sp0['smiles']}")
                    print(f"{sp1['name']}: {sp1['smiles']}")
                    print(sp0["composition"], sp1["composition"])
                    choice = int(input("Enter 0 or 1 for the first or second species: ").strip())
                    if not choice:
                        sp0["smiles"] = input(f"Enter new smiles for {sp0['name']}: ").strip()
                    else:
                        sp1["smiles"] = input(f"Enter new smiles for {sp1['name']}: ").strip()
                    print()

                if sp0.get("inchi","ZERO") == sp1.get("inchi","ONE"):
                    print(f"{sp0['name']}: {sp0['inchi']}")
                    print(f"{sp1['name']}: {sp1['inchi']}")
                    print(sp0["composition"], sp1["composition"])
                    choice = int(input("Enter 0 or 1 for the first or second species: ").strip())
                    if not choice:
                        sp0["inchi"] = input(f"Enter new inchi for {sp0['name']}: ").strip()
                    else:
                        sp1["inchi"] = input(f"Enter new inchi for {sp1['name']}: ").strip()
                    print()
    with open(f"checked-{filename}", "w") as f:
        yaml.dump(mech_data, f)

def update_atmospheric_mechanism_species(atmos_mech):
    yaml = ruamel.yaml.YAML()
    yaml.default_flow_style=False
    with open(atmos_mech, "r") as f:
        atmos_data = yaml.load(f)
    species = atmos_data["species"]
    nsp = []
    sp = species.pop(0)
    while species:
        try:
            print(sp["smiles"], sp["inchi"])
            cmps = pcp.get_compounds(sp["inchi"], "inchi")
            print(cmp.molecular_formula)
            print(cmp.isomeric_smiles)
            print(cmp.inchi)
            print(cmp.iupac_name)
            print(cmp.synonyms)
            input()
            # pop a new species only after getting through
            sp = species.pop(0)
        except Exception as e:
            print(e)
            nextsp = input("Try the next species? ")
            if nextsp:
                nsp.append(sp)
                sp = species.pop(0)

